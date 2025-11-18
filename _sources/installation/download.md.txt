## How to get the code

````{important}
For clarity, throughout this documentation variable ```OPENDIS_DIR``` is used to designate the directory that contains your downloaded OpenDiS code. For instance, you may set this directory via an environment variable ```OPENDIS_DIR``` using
```bash
export OPENDIS_DIR=${HOME}/Codes/OpenDiS.git
```
You may specify a different path that you prefer. If you use ```bash```, you can include the above line in your ```~/.bash_profile``` file so that this line is automatically executed every time you log in.
````

### Download OpenDiS together with submodules

Use the following steps to download OpenDiS into the location specified by your ```${OPENDIS_DIR}``` variable.

```bash
mkdir -p ${OPENDIS_DIR}
git clone --recurse-submodule https://github.com/OpenDiS/OpenDiS.git ${OPENDIS_DIR}
cd ${OPENDIS_DIR}
```


### Download OpenDiS and obtain submodules separately

Alternatively, you can use the following commands to achieve the same result.
```bash
mkdir -p ${OPENDIS_DIR}
git clone https://github.com/OpenDiS/OpenDiS.git ${OPENDIS_DIR}
cd ${OPENDIS_DIR}
git submodule update --init --recursive
```

### Update OpenDiS and submodules after download

If you have already downloaded OpenDiS, you can use the following commands to update it to the latest release version.
```bash
cd ${OPENDIS_DIR}
git pull
git submodule update --init --recursive
```
