# README

Author: Yifan Zhang

This implementation is based on the paper --- cuHinesBatch: Solving Multiple systems on GPUs (Ivan Martinez-Perez 
and Pedro Valero-Lara) with CUDA. 

For comparison, we also implemented serial and OpenMP version.

See Makefile for detailed information.

```sh
make all
```

After that, two executable files would be generated: `serial`, `parallel`

If any problem with the runtime, try: 
```sh
./fix-env.sh
```