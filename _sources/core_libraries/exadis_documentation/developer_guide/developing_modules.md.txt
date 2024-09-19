
## Developing modules

### General implementation approach

ExaDiS modules are typically implemented as classes which follow a generic template, e.g.:

```cpp
class MyExaDisModule : public BaseClass {
public:
    MyExaDisModule(System* system) : BaseClass(system) {
        // Initialize
    }
    void function(System* system, ...) {
        // Do some stuff
    }
    ...
};
```

where `System` is the main class that contains all information about the simulated dislocation system, including the parameters, the crystal instance, and the dislocation network object, see [System class](code_architecture.md#system-class).

#### Accessing the dislocation network

The `System` object provides the `get_serial_network()` and `get_device_network()` methods to request various instances of the dislocation network to be used in various execution spaces. These methods are implemented in the [DisNetManager class](code_architecture.md#disnetmanager-class).

Method `get_serial_network()` returns a pointer to a [SerialDisNet](code_architecture.md#serialdisnet-class) instance of the network to be used for execution in the host execution space and whose data is allocated on the host memory space.
```cpp
// Request a `SerialDisNet` instance of the dislocation network
SerialDisNet* net = system->get_serial_network();
// Perform operations in the serial execution space
...
```

Method `get_device_network()` returns a pointer to a [DeviceDisNet](code_architecture.md#devicedisnet-class) instance of the network to be used for execution in the device execution space and whose data is allocated on the device memory space or accessible from the device (e.g. via the use of unified memory).
```cpp
// Request a `DeviceDisNet` instance of the dislocation network
DeviceDisNet* net = system->get_device_network();
// Perform operation in the device execution space
...
```

```{note}
The device memory/execution spaces refer to the highest spaces available for the current build. If the code is compiled for GPU (e.g. using build option `-DKokkos_ENABLE_CUDA=On`), then the device execution/memory spaces refer to the GPU device spaces. If the code is compiled for CPU execution only, then the device  execution/memory spaces will default to the host execution/memory spaces, and any part of the code operating on a `DeviceDisNet` instance will be normally executed on the host CPU.
```

The `DisNetManager` object keeps track of the last requested instance of the dislocation network by marking it as active. When calling the `get_*_network()` methods, a memory copy is then only triggered if the requested network type instance is not currently marked as active to eliminate unnecessary memory transfers between host and device memory spaces.

For instance, in the following code example
```cpp
void function1(System* system) {
    DeviceDisNet* net = system->get_device_network();
    // do some stuff
}
void function2(System* system) {
    DeviceDisNet* net = system->get_device_network();
    // do some stuff
}
// Main body
function1(system);
function2(system);
```
no memory copy is triggered in `function2` because the `DeviceDisNet` instance was already marked as active when executing `function1`.

Conversely, in the following code example
```cpp
void function1(System* system) {
    SerialDisNet* net = system->get_serial_network();
    // do some stuff
}
void function2(System* system) {
    DeviceDisNet* net = system->get_device_network();
    // do some stuff
}
// Main body
function1(system);
function2(system);
```
a copy of the network instance from the host space to the device space will be triggered in `function2` to ensure that an up-to-date instance of the network (e.g. which may have been modified in `function1`) is accessible from the device.


### Example 1: computing the total dislocation density

As a concrete example, let's implement a function to compute the total dislocation density in the system. Calculating the dislocation density is a good starting example because it involves looping over each segment of the network, accessing their end-nodes to compute the segment length, summing up their lengths, and dividing the total length by the simulation cell volume (using appropriate units).

#### Serial implementation
As a prototype, we may start to implement this function for serial/host execution for simplicity. We could implement the following function:

```cpp
double compute_dislocation_density(System* system) {
    // Request a `SerialDisNet` instance of the dislocation network
    SerialDisNet* net = system->get_serial_network();
    
    // Loop over the local segments, compute and sum their lengths
    double density = 0.0;
    for (int i = 0; i < net->Nsegs_local(); i++) {
        auto nodes = net->get_nodes(); // generic node accessor
        auto segs = net->get_segs(); // generic segment accessor
        auto cell = net->cell;
        
        // Get segment end nodes indices
        int n1 = segs[i].n1; // end-node 1 of segment i
        int n2 = segs[i].n2; // end-node 2 of segment i
        
        // Compute segment length
        Vec3 r1 = nodes[n1].pos;
        Vec3 r2 = cell.pbc_position(r1, nodes[n2].pos); // account for PBC
        double Lseg = (r2-r1).norm();
        
        // Increment density value by segment length
        density += Lseg;
    }
        
    // Normalize by the volume and convert to 1/m^2 units
    double burgmag = system->params.burgmag; // Burgers vector magnitude in m
    density *= 1.0/net->cell.volume()/burgmag/burgmag; // 1/m^2 units
        
    // Return the dislocation density value
    return density;
}
```

#### Parallel implementation
For high-performance, we may need to implement the same function for device spaces (e.g. GPU execution). To achieve this, what we need to do is to replace the serial `for` loop by a parallel loop and operate on a `DeviceDisNet` instance of the network. Here, since we are summing the segments length into a single variable `density`, we also need to be careful to avoid race conditions (several parallel threads attempting to modify the value of the same variable will result in incorrect results). For this, Kokkos provides the `parallel_reduce()` algorithm. Our parallel implementation could look like:

```cpp
double compute_dislocation_density(System* system) {
    // Request a `DeviceDisNet` instance of the dislocation network
    DeviceDisNet* net = system->get_device_network();
    
    // Loop over the local segments, compute and sum their lengths
    double density = 0.0;
    Kokkos::parallel_reduce(net->Nsegs_local, KOKKOS_LAMBDA(const int& i, double& density_sum) {
        auto nodes = net->get_nodes(); // generic node accessor
        auto segs = net->get_segs(); // generic segment accessor
        auto cell = net->cell;
        
        // Get segment end nodes indices
        int n1 = segs[i].n1; // end-node 1 of segment i
        int n2 = segs[i].n2; // end-node 2 of segment i
        
        // Compute segment length
        Vec3 r1 = nodes[n1].pos;
        Vec3 r2 = cell.pbc_position(r1, nodes[n2].pos); // account for PBC
        double Lseg = (r2-r1).norm();
        
        // Increment density value by segment length
        density_sum += Lseg;
    }, density);
        
    // Normalize by the volume and convert to 1/m^2 units
    double burgmag = system->params.burgmag; // Burgers vector magnitude in m
    density *= 1.0/net->cell.volume()/burgmag/burgmag; // 1/m^2 units
        
    // Return the dislocation density value
    return density;
}
```
which only requires a few modifications of the serial version above.

#### Templated implementation
As a more advanced implementation, it may be beneficial to have our function working indifferently with both host/device spaces and corresponding instance types of the network. For this, we could implement the following templated version of our function:

```cpp
template<class N>
double compute_dislocation_density(System* system, N* net) {
    // Define execution policy based on network instance type
    using policy = Kokkos::RangePolicy<typename N::ExecutionSpace>;
    
    // Loop over the local segments, compute and sum their lengths
    double density = 0.0;
    Kokkos::parallel_reduce(policy(0, net->Nsegs_local), KOKKOS_LAMBDA(const int& i, double& density_sum) {
        auto nodes = net->get_nodes(); // generic node accessor
        auto segs = net->get_segs(); // generic segment accessor
        auto cell = net->cell;
        
        // Get segment end nodes indices
        int n1 = segs[i].n1; // end-node 1 of segment i
        int n2 = segs[i].n2; // end-node 2 of segment i
        
        // Compute segment length
        Vec3 r1 = nodes[n1].pos;
        Vec3 r2 = cell.pbc_position(r1, nodes[n2].pos); // account for PBC
        double Lseg = (r2-r1).norm();
        
        // Increment density value by segment length
        density_sum += Lseg;
    }, density);
        
    // Normalize by the volume and convert to 1/m^2 units
    double burgmag = system->params.burgmag; // Burgers vector magnitude in m
    density *= 1.0/net->cell.volume()/burgmag/burgmag; // 1/m^2 units
        
    // Return the dislocation density value
    return density;
}
```
which would equally work when called for serial/host execution:
```cpp
SerialDisNet* net = system->get_serial_network();
double density = compute_dislocation_density(system, net);
```
or parallel/device execution:
```cpp
DeviceDisNet* net = system->get_device_network();
double density = compute_dislocation_density(system, net);
```



### Example 2: implementing a simple force module

As another example, let's assume we want to create a new force module `ForceConstant` that sets all nodal forces to a constant vector `fval`.

#### Force base class

For force calculation modules, the base class is the `Force` class defined in `src/force.h`, and we can create the following new class:

```cpp
class ForceConstant : public Force {
private:
    // Member properties
    Vec3 fval;
    
public:
    // Force parameters
    struct Params {
        Vec3 fval = 0.0;
        Params() {}
        Params(Vec3 val) : fval(val) {}
    };
    
    // Constructor
    ForceConstant(System* system, Params params) {
        // Initialize
        fval = params.fval;
    }
    
    // Pre-compute operation
    void pre_compute(System* system) {} // nothing to pre-compute here
    
    // Global compute
    void compute(System* system, bool zero=true) {
        // Set nodal forces to fval
        ...
    }
    
    // Individual compute
    Vec3 node_force(System* system, const int& i) {
        // Compute individual force at node i
        Vec3 fi = ...
        return fi;
    }
    
    // Team individual node force implementation (optional)
    template<class N>
    KOKKOS_INLINE_FUNCTION
    Vec3 node_force(System* system, N* net, const int& i, const team_handle& team)
    {
        // Compute individual force at node i using a team of workers
        Vec3 fi = ...
        return fi;
    }
};
```

The guideline in ExaDiS is for each module to possess its own `Param` struct that serves to define the dedicated parameters of the module, and provide it to the constructor. After that, we need to implemented the various methods associated with the `Force` base class.

For a `Force` module, the 3 required methods to implement are:
- `pre_compute()`: method to perform any pre-computation that may be required to compute nodal forces. In the traditional cycle, the pre-compute function is called once at the beginning of each simulation time step. For instance, for the `ForceFFT` module, the pre-computation step is used to compute and tabulate the long-range stress field on a grid. In the `ForceSegSegList` module, the pre-computation step is used to build the list of segment pair interactions.
- `compute()`: method to perform a global computation of nodal forces, i.e. to compute the forces on all nodes of the network.
- `node_force()`: method to perform an individual force computation, i.e. to compute the force on a given node of the network. Optionally, this function can also be implemented when using a team of workers for parallel execution on device (GPU), where the individual node force computation can itself be parallelized across several concurrent threads using hierarchical parallelization.

#### Serial implementation

In our simple example, there is no need to perform any pre-computation, so we can leave the `pre_compute()` method empty. The next step is to implement the `compute()` method. Recall we want to set the nodal forces to a constant value. As a first attempt/prototype, we may want to implement this method in serial on the CPU for simplicity. Here, we could do something like this:

```cpp
// Global compute (serial implementation)
void compute(System* system, bool zero=true) {
    // Request a `SerialDisNet` instance of the dislocation network
    SerialDisNet* net = system->get_serial_network();
    
    // Zero out the nodal forces if requested (mandatory instruction)
    if (zero) zero_force(net);
    
    // Set nodal forces to fval using a simple for loop
    for (int i = 0; i < net->Nnodes_local; i++) {
        auto nodes = net->get_nodes(); // generic node accessor
        nodes[i].f = fval;
    }
}
```

where we loop over the local nodes and assign `fval` to their force property. We're done, and we can check with a test case that this is working properly. 

#### Parallel implementation

Now, let's imagine that we are happy with our implementation, but want to use this module in a production run on GPU. In this case, executing the `compute()` method in the serial execution space is going to be very inefficient. (This would be the case even for such a trivial force module. This is because other force/mobility modules used in the simulation will likely be executed on GPU, hence the call to the `system->get_serial_network()` will trigger memory copies from the GPU to the CPU, while the `system->get_device_network()` in the other GPU modules will trigger the reverse memory copies. These back-and-forth movements may significantly hit the performance.) Alternatively, we can now implement a new version of the `compute()` method that will be executed on the device space. Here, the implementation is trivial because each index `i` of the loop is independent and thus can be parallelized. As such, all we need to do is to replace the serial `for` loop with a `Kokkos::parallel_for` loop, while now operating on a `DeviceDisNet` object:

```cpp
// Global compute (parallel implementation)
void compute(System* system, bool zero=true) {
    // Request a `DeviceDisNet` instance of the dislocation network
    DeviceDisNet* net = system->get_device_network();
    
    // Zero out the nodal forces if requested (mandatory instruction)
    if (zero) zero_force(net);
    
    // Set nodal forces to fval using a parallel_for loop
    Vec3 f = fval; // set vector variable in local scope
    Kokkos::parallel_for(net->Nnodes_local, KOKKOS_LAMBDA(const int i) {
        auto nodes = net->get_nodes();
        nodes[i].f = f;
    });
    Kokkos::fence();
}
```

Here we are using a lambda expression to define the parallel kernel using the `KOKKOS_LAMBDA` macro. Note that the vector `fval` is copied to a local vector `f` in the method scope in order to avoid implicit capture of a class member in the lambda function. As slightly more advanced implementations, we could use class `operator()` or a dedicated functor to define the same parallel kernel without using lambda functions. Note that we also need to call the `Kokkos::fence()` function to make sure we have finished executing the potentially asynchronous `parallel_for` kernel launch before leaving the method.

Finally, we can also implement the `node_force()` method to compute the force on a single node. Here this is trivial and we can just have:

```cpp
// Individual compute
Vec3 node_force(System* system, const int& i) {
    // Compute individual force at node i
    Vec3 fi = fval;
    return fi;
}
```

and we can implement the same for the `node_force()` method using team workers. The team implementation is optional but allows other modules that use `node_force()` calls in parallel kernels (e.g. module `TopologyParallel`) to be used with our `ForceConstant` module.


### Example 3: implementing a segment-based force module

As a slightly more complex example, let's now imagine that we want to implement a force module `ForceSegConstant` for which the force on each segment is the constant `fval` multiplied by the length of the segment, and that the segment force is equally distributed between the two end-nodes.

#### Serial implementation

In a serial fashion, we could write the `compute()` function as:

```cpp
// Global compute (serial implementation)
void compute(System* system, bool zero=true) {
    // Request a `SerialDisNet` instance of the dislocation network
    SerialDisNet* net = system->get_serial_network();
    
    // Zero out the nodal forces if requested (mandatory instruction)
    if (zero) zero_force(net);
    
    // Loop over segments to aggregate forces at the nodes
    for (int i = 0; i < net->Nsegs_local; i++) {
        auto nodes = net->get_nodes(); // generic node accessor
        auto segs = net->get_segs(); // generic segment accessor
        auto cell = net->cell;
        
        // Get segment end nodes indices
        int n1 = segs[i].n1; // end-node 1 of segment i
        int n2 = segs[i].n2; // end-node 2 of segment i
        
        // Compute segment length
        Vec3 r1 = nodes[n1].pos;
        Vec3 r2 = cell.pbc_position(r1, nodes[n2].pos); // account for PBC
        double L = (r2-r1).norm();
        
        // Distribute segment force at nodes
        nodes[n1].f += 0.5*L*fval;
        nodes[n2].f += 0.5*L*fval;
    }
}
```

Here, we now loop over the local segments, compute the segments length from their end-nodes positions, and distribute the segment force equally by incrementing the end-nodes force values.

#### Parallel implementation

We can notice that if we want to implement the same method in a parallel fashion, we would need to be careful. This is because indices `i` in the loop are no longer fully independent: two distinct segments `i` may need to access some of the same nodes `n1` or `n2`. When running the above loop kernel in parallel, this may create race conditions and yield incorrect results.

To simplify the implementation, ExaDiS provides additional base classes that abstract away this type of complexity. In this particular case, base class `ForceSeg` in `src/force.h` provides a simple way to implement our desired parallel kernel without having to worry about the parallelism aspect of it. All what base class `ForceSeg` requires is the kernel inside of the parallel loop to be provided in the form of a `segment_force()` method that returns end-nodes forces. To provide it, we simply need to define a struct that implements our desired `segment_force()` kernel:

```cpp
struct SegConstant
{
    // Flags to instruct what kernels are implemented in the struct
    static const bool has_pre_compute = false;
    static const bool has_compute_team = false;
    static const bool has_node_force = false;
    
    Vec3 fval;
    
    // Force parameters
    struct Params {
        Vec3 fval = 0.0;
        Params() {}
        Params(Vec3 val) : fval(val) {}
    };
    
    // Constructor
    SegConstant(System* system, Params params) {
        // Initialize
        fval = params.fval;
    }
    
    // Segment force kernel
    template<class N>
    KOKKOS_INLINE_FUNCTION
    SegForce segment_force(System* system, N* net, const int& i) 
    {
        auto nodes = net->get_nodes(); // generic node accessor
        auto segs = net->get_segs(); // generic segment accessor
        auto cell = net->cell;
        
        // Get segment end nodes indices
        int n1 = segs[i].n1; // end-node 1 of segment i
        int n2 = segs[i].n2; // end-node 2 of segment i
        
        // Compute segment length
        Vec3 r1 = nodes[n1].pos;
        Vec3 r2 = cell.pbc_position(r1, nodes[n2].pos); // account for PBC
        double L = (r2-r1).norm();
        
        Vec3 f1 = 0.5*L*fval;
        Vec3 f2 = 0.5*L*fval;
        return SegForce(f1, f2);
    }
    
    static constexpr const char* name = "SegConstant";
};
```

and then declare our `ForceSegConstant` force module as the base class `ForceSeg` templated with our struct `SegConstant` implementing the kernel:

```cpp
typedef ForceSeg<SegConstant> ForceSegConstant;
```

In struct `SegConstant`, all we had to do is to pretty much copy the inside of our serial loop into the `segment_force()` method. We have also templated the network instance type `N` so that the same kernel can be compiled for serial or device execution spaces using indifferently `SerialDisNet` or `DeviceDisNet` instances of the network. When compiled for GPU, the segment forces will be computed in a highly parallel fashion on the device (GPU) space by default, and forces at nodes will be aggregated properly avoiding race conditions, following the machinery implemented in base class `ForceSeg` and here abstracted from the user. In addition, when using this approach the associated `node_force()` method becomes automatically created as well, without us having to explicitly define it. However, if we want to implement a dedicated `node_force()` method, we could also do that by setting flag `has_node_force = true` and implementing method `node_force()` in our `SegConstant` struct, in which case its base class implementation will be overridden. Similarly, we could provide a team implementation or pre-compute methods.

In the code, base class `ForceSeg` is for instance used to compute the core force in `src/force_types/force_lt.h`, or implement the N^2 force model in `src/force_types/force_n2.h`.
