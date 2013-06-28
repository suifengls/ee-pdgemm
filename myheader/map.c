#include "map.h"
#include <stdio.h>
#include <stdlib.h>

// in binding=rr mode
void mapping(int mycol, int type)
{
   char H_Freq[7]={"2500000"};
//   char M_Freq[7]={"1800000"};
   char L_Freq[6]={"800000"};
   FILE *fp;
   char *DVFS=L_Freq;

   if(type==1) DVFS=H_Freq;
   else if(type==0) DVFS=L_Freq;

   if(mycol==0)
   {
      fp = fopen("/sys/devices/system/cpu/cpu0/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==1)
   {
      fp = fopen("/sys/devices/system/cpu/cpu1/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==2)
   {
      fp = fopen("/sys/devices/system/cpu/cpu2/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==3)
   {
      fp = fopen("/sys/devices/system/cpu/cpu3/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==4)
   {
      fp = fopen("/sys/devices/system/cpu/cpu4/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==5)
   {
      fp = fopen("/sys/devices/system/cpu/cpu5/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==6)
   {
      fp = fopen("/sys/devices/system/cpu/cpu6/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==7)
   {
      fp = fopen("/sys/devices/system/cpu/cpu7/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
}


void mapping2(int myrow, int mycol, int type)
{
   char H_Freq[7]={"2500000"};
//   char M_Freq[7]={"1800000"};
   char L_Freq[6]={"800000"};
   FILE *fp;
   char *DVFS=L_Freq;

   if(type==1) DVFS=H_Freq;
   else if(type==0) DVFS=L_Freq;

if(myrow < 4)
{ 
   mycol=mycol%4;

   if(mycol==0)
   {
      fp = fopen("/sys/devices/system/cpu/cpu0/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==1)
   {
      fp = fopen("/sys/devices/system/cpu/cpu1/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==2)
   {
      fp = fopen("/sys/devices/system/cpu/cpu2/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==3)
   {
      fp = fopen("/sys/devices/system/cpu/cpu3/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
}
else
{ 
   if(mycol < 4)
      mycol = mycol + 4;
   
   if(mycol==4)
   {
      fp = fopen("/sys/devices/system/cpu/cpu4/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==5)
   {
      fp = fopen("/sys/devices/system/cpu/cpu5/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==6)
   {
      fp = fopen("/sys/devices/system/cpu/cpu6/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
   else if(mycol==7)
   {
      fp = fopen("/sys/devices/system/cpu/cpu7/cpufreq/scaling_setspeed","w");
      fwrite(DVFS,6+type,sizeof(char),fp);
      fclose(fp);
   }
}

}
