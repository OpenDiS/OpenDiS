## How to get the code

### Download OpenDiS together with submodules

For clarity, let us set the environment variable ```OPENDIS_DIR``` to be the directory that will contain your downloaded OpenDiS code.  For example, you may use
```bash
export OPENDIS_DIR=${HOME}/Codes/OpenDiS.git
```
But you may specify a different path that you prefer.  If you use ```bash```, you can include the above line in your ```~/.bash_profile``` file so that this line is automatically executed every time you log in.

Use the following steps to download OpenDiS into the location specified by your ```${OPENDIS_DIR}``` variable.

```bash
mkdir -p ${OPENDIS_DIR}
git clone --recurse-submodule https://gitlab.com/micronano/OpenDiS.git ${OPENDIS_DIR}
cd ${OPENDIS_DIR}
```


### Download OpenDiS and obtain submodules separately

Alternatively, you can use the following commands to achieve the same result.
```bash
mkdir -p ${OPENDIS_DIR}
git clone https://gitlab.com/micronano/OpenDiS.git ${OPENDIS_DIR}
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
