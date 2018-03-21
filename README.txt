Name: Weiyu Zhu
UT EID: wz4245

To build, use the following commands:

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8

excutable file will be "build/bin/menger"

The saved file "geometry.obj" should be saved in the folder where the executable file is executed, e.g. "bin".

All the required parts are completed.

For the optional parts:

1. Constant-width wireframes are enerated.

2. In the current position of the light source, I draw a sphere there, however, it's just white, not glowing =_=

3. I create a "skybox" that wraps my scene.   To show the skybox, the executable file "menger" should be run from folder "bin", the image files are in directory "src/cubemap".