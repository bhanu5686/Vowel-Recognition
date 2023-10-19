#include "stdafx.h"   

#define ld long double 
ld raisedSineWeight[12]={0.0};
char vowel[6]={'a','e','i','o','u'};
ld tokhuraWeight[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
ld PI =  3.141592653589793238; 


bool ReadFile(ld *samples,const char  *fileName)
{
    FILE *fp;
    fp = fopen(fileName, "r");
    if(!fp)
        return false;
    else
    {
        int i=5,a;
        char buffer[100];
        while(i--)
        {
            fgets(buffer, 100, fp);
        }
        i=0;
        while(fscanf(fp,"%d",&a)!=EOF)
        {
            samples[i++]=a;
        }        
    }
    fclose(fp);
    return true;
}


bool WriteFile(ld avgCIS[5][13],const char  *ReffileName)
{
    FILE *fp;
    fp = fopen(ReffileName, "w");
    if(!fp)
    {
        return false;
    }
    else
    {
        int i,j;
        for(i=0; i<5 ;i++)
        {
            for(j=1;j<13 ;j++)
            {
                fprintf(fp, "%Lf  ", avgCIS[i][j]);
            }
            fprintf(fp, "\n");
        }               
    }

    fclose(fp);
    return true;
}


bool ReadStablePart(ld *stableSamples, ld *samples,int size)
{    
    ld STE = 0.0;   
    int end = 320*5;
  
    if(size<end)
    {        
        return false;
    }
  
    int start = 0;
  
    for(int i=0;i<end;i++)
    {      
        STE += (samples[i]*samples[i]);
    }
  
    for(int i=end;i+end<size;i+=end)
    {    
        
        ld currSTE = 0.0;
    

        for(int j=i;j<i+end;j++)
        {
            currSTE += (samples[j]*samples[j]);
        }
    
        if(currSTE > STE)
        {
            STE = currSTE;
            start = i;
        }
    }
  
    int x=0;
    for(int i=start;i<start+end;i++)
    {           
        stableSamples[x++] = samples[i];
    }
    return true;
}


void DCShift(ld *samples,int size)
{    
    int i=0;
    ld sum=0;
    for(i=0;i<size;i++)
    {
        sum = sum + samples[i];
    }
    sum = sum/(ld)size ;
    for(i=0;i<size;i++)
    {
        samples[i] = samples[i] - sum ;
    }    
    return;
}


void Normalization(ld *samples,int size)
{

    ld maxx = 0.0;
    int i=0;  

    for(;i<size;i++)
    {
        if(maxx < samples[i])
            maxx = samples[i];  
    }  

    for(int i=0;i<size;i++)
    {
        samples[i] /= maxx;        
        samples[i] *= 5000;
    }
    return;
}


void ComputeRIS(ld *stableSamples,ld *RIS, int start, int end)
{
    int i,j,x=320;
    ld n = (ld)x;
    for(i=0;i<=12;i++)  
    {
        ld sum = 0.0;
        for(j=start;j+i<end;j++)
        {
            sum = sum + (stableSamples[j] * stableSamples[i+j]);
        }
        sum = sum/n;
        
        RIS[i] = sum;
    }
    return;
}



void ComputeAIS(ld *RIS,ld *AIS)
{
    ld energy[13] = {0.0};
    ld K[13] = {0.0};
    ld newAlpha[13] = {0.0};
    ld prevAlpha[13] = {0.0};

    energy[0] = RIS[0];
    
    int i=0,j=0;
   
    for(i=1;i<=12;i++)
    {
        if(i==1)
            K[i] = RIS[i]/RIS[0]; 
        else 
        {
            ld sum=0.0;
            for(j=1;j<i;j++)
            {
                prevAlpha[j] = newAlpha[j]; 
                sum += (prevAlpha[j] * RIS[i-j]);
            }
            K[i] = (RIS[i] - sum)/energy[i-1];
        }

        newAlpha[i] = K[i]; 
        for(j=1;j<i;j++)
        {
            newAlpha[j] = prevAlpha[j] - (K[i] * prevAlpha[i-j]); 
        } 
        energy[i] = (1 - (K[i] * K[i])) * energy[i-1]; 
    }
    for(i=1;i<=12;i++)
        AIS[i] = newAlpha[i];
    return;
}


void ComputeCIS(ld *AIS, ld *RIS, ld *CIS)
{
    CIS[0] = logl(RIS[0]);
    
    int m=0,k=0;
    for(m=1;m<=12;m++)
    {        
        ld sum = 0.0;
        
        for(k=1;k<m;k++)
        {
            sum += ((((ld)k)*CIS[k]*AIS[m-k])/((ld)m)); 
        }
        
        ld value = AIS[m] + sum;
        CIS[m] = value;
    }
    return;
}


void RiasedSineWindow(ld *CIS)
{
    for(int i=1;i<=12;i++)
    {
        CIS[i] = raisedSineWeight[i-1] * CIS[i];    
    }  
  return;
}



void ComputeRSW(ld *raisedSineWeight)
{
    for(int i=1;i<=12;i++)
    {    
        ld m = (ld)i;    
        ld theta = (PI * m)/12.0;
        ld value = sin(theta);
        value = 1.0 + (6.0 * value);
        raisedSineWeight[i-1] = value;
    }
    return;
}


void ComputeAvgAllCIS(ld allCIS[50][13], ld avgCIS[5][13]) 
{
    int i,j;
    for(i=0;i<(5*10);i++)
    {    
        for(j=1;j<=12;j++)
        {
            int x = i%5;
            int y = j;
            avgCIS[x][y] += allCIS[i][j]; 
        }
    }    
        
    for(i=0;i<5;i++)
    {
        for(j=1;j<=12;j++)
        {
            avgCIS[i][j] /= 10.0;
        }
    }

    return;
}

ld TokhuraDistance(ld allCIS[50][13], ld C[5][13])
{
  
  ld tDistance = LDBL_MAX;
  
  for(int i=0;i<5;i++)
  {
    
    ld distance = 0.0;
    
    for(int j=1;j<13;j++)
    {
      ld difference = allCIS[i][j] - C[i][j];       
      distance += (tokhuraWeight[j-1]*(difference * difference));
    }
    
    
    if(distance < tDistance)
    {
      tDistance = distance;
    }
  }
  return tDistance;
}


char VowelRecognize(ld allCIS[50][13])
{

    char output = '.';
    ld minTDistance = LDBL_MAX;
    char vowels[5] = {'a', 'e', 'i', 'o', 'u'};
    ld C[5][13] = {0.0}; 
    string fileName = "";
  
    for(int i=0;i<5;i++)
    {         
        string s = "";
        s+=vowels[i];
        string ReffileName = ".\\reference_files\\" + s + ".txt";

        FILE *fp;
        fp = fopen(ReffileName.c_str(), "r");
        if(!fp)
            return false; 
        else
        {    
            int j,k;
            ld value = 0.0;
            for(j=0;j<5;j++)
            {                  
                for(k=1;k<=12;k++)
                {
                    fscanf(fp,"%Lf",&value);
                    C[j][k] = value;
                }
            }
            ld tDistance = TokhuraDistance(allCIS, C);
            if(tDistance < minTDistance)
            {  
                minTDistance = tDistance;
                output = vowels[i];
            }
    
            fclose(fp);            
        }  
    }  
    return output;
}


int sizeOfFile(const char  *fileName)
{
    FILE *fp;
    int cnt=0;
    fp = fopen(fileName, "r");
    if(!fp)
        return -1;
    else
    {
        int i=5,a;
        char buffer[100];
        while(i--)
        {
            fgets(buffer, 100, fp);
        }
        while(fscanf(fp,"%d",&a)!=EOF)
        {
            cnt++;
        }        
    }
    fclose(fp);
    return cnt;
}


int _tmain(int argc, _TCHAR* argv[])
{
    ComputeRSW(raisedSineWeight);

    int i=0,j=0,k=0,x=0;    
    printf("Reference files generation starts..\n\n");
    for(i=0;i<5;i++)
    {
        x=0;      

        ld avgCIS[5][13] = {0.0}; 
        ld allCIS[50][13] = {0.0};  

        for(j=1;j<=10;j++)
        {            
            string s = "";
            s+=vowel[i];
            string fileName = ".\\input_files\\" + s + "_" + to_string((long long int)j) + ".txt";
            int no = sizeOfFile(fileName.c_str());
            if(no == -1)
            {
                printf("unable to open file");
                return 0;
            }
            ld* samples;
            samples = (ld*)malloc(no * sizeof(ld));

            if(!ReadFile(samples,fileName.c_str()))
            {
                printf("unable to open input_files!..");
                return 0;
            }
            DCShift(samples,no);
            Normalization(samples,no);

            ld* stableSamples;
            stableSamples = (ld*)malloc(1600 * sizeof(ld));

            if(!ReadStablePart(stableSamples, samples ,no))
            {                
                printf("Recording is too small.");
                return 0;
            }


            for(int k=0;k<5;k++)
            {   
                ld RIS[13]={0.0};
                ld AIS[13]={0.0};
                ld CIS[13]={0.0};
                int start = (320 * k);
                int end = (320 * (k+1));
                ComputeRIS(stableSamples,RIS, start, end);          
                
                ComputeAIS(RIS,AIS);   
                
                ComputeCIS(AIS,RIS,CIS);    
                
                RiasedSineWindow(CIS);  
                 
                for(int n=0;n<13;n++)
                {
                    allCIS[x][n] = CIS[n];
                }
                x++;
            }
			free(samples);
			free(stableSamples);
        }
        
        ComputeAvgAllCIS(allCIS,avgCIS);      
        

        string st = "";
        st+=vowel[i];
        string ReffileName = ".\\reference_files\\" + st + ".txt";
   
        if(!WriteFile(avgCIS, ReffileName.c_str()))
        {
            printf("unable to open Reference_file!..");
            return 0;
        }        
        
    }

    printf("Reference files is generated.\n\n");
    printf("Testing starts...\n\n");

    x=0;
    for(i=0;i<5;i++)
    {   
            
        for(j=11;j<=20;j++)
        {
            string s = "";
            s+=vowel[i];
            string fileName = ".\\testing_files\\" + s + "_" + to_string((long long int)j) + ".txt";            
           
            ld allCIS[5][13] = {0.0}; 
            x=0; 
            int no = sizeOfFile(fileName.c_str());
            if(no == -1)
            {
                printf("unable to open file");
                return 0;
            }
            ld* samples;
            samples = (ld*)malloc(no * sizeof(ld));

            if(!ReadFile(samples, fileName.c_str()))
            {      
                printf("unable to open input_file ");
                return 0;
            }

            DCShift(samples,no);   
            Normalization(samples,no);

            ld* stableSamples;
            stableSamples = (ld*)malloc(1600 * sizeof(ld));

            if(!ReadStablePart(stableSamples, samples, no))
            {     
                printf("Recording is too small");
                return 0;
            }    
      
            for(int k=0;k<5;k++)
            {   
                ld RIS[13]={0.0};
                ld AIS[13]={0.0};
                ld CIS[13]={0.0};
                int start = (320 * k);
                int end = (320 * (k+1));
                ComputeRIS(stableSamples,RIS, start, end);
                ComputeAIS(RIS,AIS);
                ComputeCIS(AIS,RIS,CIS);
                RiasedSineWindow(CIS);     
                for(int n=0;n<13;n++)
                {
                    allCIS[x][n] = CIS[n];
                }
                x++;
            }

            char vow = VowelRecognize(allCIS);    
            if(vow != '.')
            {
                cout<<"Vowel present in file - "<<fileName<<" is - "<<vow<<endl;
            }
			free(samples);
			free(stableSamples);
        }
        printf("\n");
    }    
    return 0;
}

