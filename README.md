Energy Efficient - Parallel Matrix-Matrix Multiplication
========
PS: You should have the permission to change the cpu frequency by writing /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed
Files for this project
--------
1. demon

   contains the main function files
2. demon.pbs

   PBS file to hand in jobs
3. hostfile
 
   hostfile for MPI
4. myheader
 
   contains pipeline, mypdgemm, mapping functions
5. sample_output.dat
 
   a output file of eepdgemm

---
Paper: [Energy Efficient Parallel Matrix-Matrix Multiplication for DVFS-Enabled Clusters](http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=6337486&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D6337486)
