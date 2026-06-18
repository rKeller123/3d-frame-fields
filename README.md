# 3d-frame-fields

## building in linux

```bash
# release mode

cd build
cmake ..
cmake --build .

./3d-odec-frame-fields

# debug
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
cmake --build .

flamegraph -- ./3d-odeco-frame-fields 10
xdg-open flamegraph.svg
```