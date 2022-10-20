#ifndef _OPENDIS_H
#define _OPENDIS_H

class OPENDIS {
  public:
    class DisNetwork *network;  // dislocation network structure

    class KokkosOds *kokkos;    // KOKKOS accelerator class
    class Python *python;       // Python interface

    OPENDIS();
    ~OPENDIS();

}

#endif

