#include <iostream>
#include <stdlib.h>
#include "cmath"
#include <ctime>
#include <fstream>
using namespace std;
class   AirShower{
	private: float EnergyParticle;
	private: float XFlowAxis;
	private: float YFlowAxis;
	private: float Teta;
	private: float Fi;
	private: float ValueInPoint;
	private: float DegreeDrop;
	private: float ValueMax;
	
	private: float Rand(float Max){
		float e;
        e=rand() % 10000;
        e=e * Max / 10000;
        return e;
	};	
	public:  float CalculatedValueShowerInPoint(float XCoord,float YCoord){
	float XDifference=(XFlowAxis - XCoord);
	float YDifference=(YFlowAxis - YCoord);	
	float	Distance = sqrt ( pow(XDifference,2) + pow(YDifference,2) );
	float ValuePoint = 10 / (ValueMax + pow(Distance/10,2.5))+ Rand(1 / (ValueMax + pow(Distance/10,2.2)));
	return ValuePoint;
	};
	AirShower(float XFlowAxisInd, float YFlowAxisInd,float Energy,float TetaInd,float FiInd){
		ValueMax=Energy;
		XFlowAxis=XFlowAxisInd;
		YFlowAxis=YFlowAxisInd;
		Teta=TetaInd;
		Fi=FiInd;
	};
};
class FieldDetectors{	
	private: int SumDetectorLength;
	private: int SumDetectorWigth;
	private: float Xmax;
	private: float Ymax;
	private: float **MatrixRec;
	private: float DistanceBetweenDetector;
 	private: float XCoordinates;
 	private: float YCoordinates;
	private: float Rand(float Max){
		float e;
        e=rand() % 10000;
        e=e * Max / 10000;
        return e;
	};	
	public: float GetXmax(){
 		return Xmax;
	 };
	public: float GetYmax(){
 		return Ymax;
	 };
	public: float SetXmax(float m){
 		Xmax=m;
	 };
	 public: float SetYmax(float m){
 		Ymax=m;
	 }; 
 	
 	public: float GetXCoordinates(){
 		return XCoordinates;
	 };
 	public: float GetYCoordinates(){
 		return YCoordinates;
	 };
	public: float GetSumDetectorLength(){
 		return SumDetectorLength;
	 };
	public: float GetSumDetectorWigth(){
 		return SumDetectorWigth;
	 };
 
	public: void SetValueInMatrix(float Value,int NLength, int NWigth ){
		MatrixRec[NLength][NWigth]=Value;
	};	

    public: void AddValueInMatrix(float Value,int NLength, int NWigth ){
		MatrixRec[NLength][NWigth] = MatrixRec[NLength][NWigth] + Value;
	};	

	public: float FoundCoordinatesRectangular(int NumberLength, int NumberWigth){
		XCoordinates=NumberLength * DistanceBetweenDetector;
		YCoordinates=NumberWigth  * DistanceBetweenDetector;	
	};
	public: float FoundCoordinatesHexagon(int NumberLength, int NumberWigth){
		XCoordinates=NumberLength * DistanceBetweenDetector;
		YCoordinates=NumberWigth * DistanceBetweenDetector;
	}
	
	public: void CreateMatrixRectangular() {
	    MatrixRec=new float*[SumDetectorLength]; 
       	    for (int NLength=0; NLength<SumDetectorLength; NLength++)  
		   		MatrixRec[NLength]=new float[SumDetectorWigth]; 
	};
	
	
	public: void GenerateNoise(float Max){
	    for(int NLength=0;NLength<SumDetectorLength;NLength++){
			for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++){
					MatrixRec[NLength][NWigth]=Rand(Max);
			}
		}			
	};
	public: void CoutMatrix(){	
		for(int NLength=0;NLength<SumDetectorLength;NLength++){
			for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++){
				cout<<MatrixRec[NLength][NWigth]<<" ";
			}
			cout<<endl;
		}  
	};
	public: void FoutMatrix(){
		ofstream DeleteFail;DeleteFail.open("C:\\Users\\1\\Documents\\данные2.txt",ios::out);DeleteFail.close();
  		ofstream L; 
  		L.open("C:\\Users\\1\\Documents\\данные2.txt",ios::app);
		for(int NLength=0;NLength<SumDetectorLength;NLength++){
			for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++){
				L<<MatrixRec[NLength][NWigth]<<"  ";
			}
			L<<endl;
		}  
	}
	public: void FindAxisMaxDetector(){
		float Maximum;
        	for(int NLength=0;NLength<SumDetectorLength;NLength++) {
				for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++) {
					if (MatrixRec[NLength][NWigth]>Maximum){
						Maximum=MatrixRec[NLength][NWigth];
						Xmax=NLength;
						Ymax=NWigth;
					}
				}
			}
	};	
	
	public: void FindAxisMNC(float F){
		float SumEX;
	    float SumEY;
	    float SumE,R;
	    SumEX=0;  SumEY=0;  SumE=0;
		for(int NLength=0;NLength<SumDetectorLength;NLength++)	 {
			for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++)    {
				FoundCoordinatesRectangular( NLength , NWigth );	
				SumEX=SumEX+pow(MatrixRec[NLength][NWigth],2) * XCoordinates;
				SumEY=SumEY+pow(MatrixRec[NLength][NWigth],2) * YCoordinates;
				SumE=SumE+pow(MatrixRec[NLength][NWigth],2);
			}
		};
		Xmax=SumEX / SumE;
		Ymax=SumEY / SumE;
		R=sqrt(pow((5-Xmax),2)+pow((5-Ymax),2));
		if (4.51-Xmax>0 )
		{
		Xmax=Xmax-(F*0.03+0.19)*pow(abs(4.51-Xmax),2.2);}
		else{
		Xmax=Xmax+(F*0.03+0.19)*pow(abs(4.51-Xmax),2.2);
	}	if (4.51-Ymax>0 )
		{
		Ymax=Ymax-(F*0.03+0.19)*pow(abs(4.51-Ymax),2.2);}
		else{
		Ymax=Ymax+(F*0.03+0.19)*pow(abs(4.51-Ymax),2.2);
	}			
	};
		public: void FindAxisMNCHex(){
		float SumEX;
	    float SumEY;
	    float SumE;
	    SumEX=0;  SumEY=0;  SumE=0;
		for(int NLength=0;NLength<SumDetectorLength;NLength++)	 {
			for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++)    {
				FoundCoordinatesHexagon( NLength , NWigth );	
				SumEX=SumEX+pow(MatrixRec[NLength][NWigth],2) * XCoordinates;
				SumEY=SumEY+pow(MatrixRec[NLength][NWigth],2) * YCoordinates;
				SumE=SumE+pow(MatrixRec[NLength][NWigth],2);
			}
		};
		Xmax=SumEX / SumE;
		Ymax=SumEY / SumE;
		Xmax=Xmax-0.55*(9.5-Xmax);
		Ymax=Ymax-0.35*(4.5-Ymax);	
	};
		public: float FindDif(float Xm, float Ym,float ValueMax,float Power){ 	 float Dif=0;
			for(int NLength=0;NLength<SumDetectorLength;NLength++)	 {
			for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++)    {
			float XCoord=NLength * DistanceBetweenDetector;
		    float YCoord=NWigth  * DistanceBetweenDetector;	
		    float XDifference=(Xm - XCoord);
			float YDifference=(Ym - YCoord);	
			float	Distance = sqrt ( pow(XDifference,2) + pow(YDifference,2) );
			float Value = 10 / (ValueMax + pow(Distance/10,Power));
			Dif=Dif + abs(Value-MatrixRec[NLength][NWigth]);
			}
		  }
		  return Dif;
		}
        public: CrossSection(float YAxis,float A){
			for(int NLength=0;NLength<SumDetectorLength;NLength++)	 {
			for(int NWigth=0;NWigth<SumDetectorWigth;NWigth++)    {
float V1,Y1,V2,S1,S2,X1,SumLeft,SumRight;
FoundCoordinatesRectangular(NLength,NWigth);
if((NLength<SumDetectorLength-1)&&(NWigth<SumDetectorWigth-1)){
V1=pow(DistanceBetweenDetector,2)*(MatrixRec[NLength][NWigth]*(1-1/sqrt(2))+MatrixRec[NLength+1][NWigth]/2/sqrt(2)+MatrixRec[NLength][NWigth+1]/2/sqrt(2));
Y1=YAxis+tan(A)*XCoordinates;
if((YCoordinates-Y1)/tan(A)<=(XCoordinates+DistanceBetweenDetector)){
S1=pow((YCoordinates-Y1),2)/(1+tan(A))/2;
SumLeft+=pow(DistanceBetweenDetector,2)/2/S1*V1;
SumRight+=pow(DistanceBetweenDetector,2)/2/(pow(DistanceBetweenDetector,2)/2-S1)*V1;
}
else{
X1=(Y1-YCoordinates+XCoordinates*(1-tan(A)))/(1-tan(A));
S2=(Y1-YCoordinates-DistanceBetweenDetector)*(sqrt(pow(X1,2)+pow((YCoordinates-DistanceBetweenDetector+X1-XCoordinates-DistanceBetweenDetector),2)))/4*sqrt(2);
SumRight+=pow(DistanceBetweenDetector,2)/2/S2*V1;
SumLeft+=pow(DistanceBetweenDetector,2)/2/(pow(DistanceBetweenDetector,2)/2-S2)*V1;
}
}    
if((NLength>SumDetectorLength+1)&&(NWigth>SumDetectorWigth+1)){
	V2=pow(DistanceBetweenDetector,2)*(MatrixRec[NLength][NWigth]*(1-1/sqrt(2))+MatrixRec[NLength-1][NWigth]/2/sqrt(2)+MatrixRec[NLength][NWigth-1]/2/sqrt(2));
	Y1=YAxis+tan(A)*XCoordinates;
	if((YCoordinates-Y1)/tan(A)<=(XCoordinates+DistanceBetweenDetector)){
	S1=pow((YCoordinates-Y1),2)/(1+tan(A))/2;
	SumLeft+=pow(DistanceBetweenDetector,2)/2/S1*V1;
	SumRight+=pow(DistanceBetweenDetector,2)/2/(pow(DistanceBetweenDetector,2)/2-S1)*V1;
	}
	else{
	X1=(Y1-YCoordinates+XCoordinates*(1-tan(A)))/(1-tan(A));
	S2=(Y1-YCoordinates-DistanceBetweenDetector)*(sqrt(pow(X1,2)+pow((YCoordinates-DistanceBetweenDetector+X1-XCoordinates-DistanceBetweenDetector),2)))/4*sqrt(2);
	SumRight+=pow(DistanceBetweenDetector,2)/2/S2*V1;
	SumLeft+=pow(DistanceBetweenDetector,2)/2/(pow(DistanceBetweenDetector,2)/2-S2)*V1;
	}
}
}
}
}
		
	FieldDetectors(int Length, int Wigth,float DistanceBetween){
		SumDetectorLength = Length;
		SumDetectorWigth =  Wigth;
		DistanceBetweenDetector=DistanceBetween;
		}	
};


int main(){ 
ofstream D; D.open("C:\\Users\\1\\Documents\\данные.txt");D.close();
	int I;
	float f=1;
  	ofstream L,L1; 
  	float Field,Xmaxsimum,Ymaxsimum,Sumu,H,M,N,POWER,g1,s1=0;
	srand(time(0));
	L.open("C:\\Users\\1\\Documents\\данные.txt",ios::app);
	L1.open("C:\\Users\\1\\Documents\\данные1.txt",ios::app);
for(float k=49;k<50;k++){
for(float i=49;i<50;i++){
	FieldDetectors firstfile(10,10,1);
			AirShower ProtonAirShower(i/10,k/10,0.1,0,0);
			firstfile.CreateMatrixRectangular();	
			firstfile.GenerateNoise(f);
			for(int NLength=0;NLength<firstfile.GetSumDetectorLength();NLength++) {
				for(int NWigth=0;NWigth<firstfile.GetSumDetectorWigth();NWigth++) {
					float Value;
					firstfile.FoundCoordinatesRectangular(NLength, NWigth);
					//Value=ProtonAirShower.CalculatedValueShowerInPoint(firstfile.GetXCoordinates(),firstfile.GetYCoordinates());
					//firstfile.AddValueInMatrix(Value,NLength,NWigth);		
				}	
			}
M=0.5;
while(M<=15){
M+=1;
firstfile.CrossSection(M,0.00001);
cout<<M<<endl;
}
	
	//s1=sqrt(pow((i/10-firstfile.GetXmax()),2)+pow((k/10-firstfile.GetYmax()),2));
	//L<<s1<<endl;
	//cout<<sqrt(pow((i/10-firstfile.GetXmax()),2)+pow((k/10-firstfile.GetYmax()),2))<<"  "<<firstfile.GetXmax()<<"  "<<firstfile.GetYmax()<<"  "<<endl;
	//L<<sqrt(pow((i/10-firstfile.GetXmax()),2)+pow((k/10-firstfile.GetXmax()),2))<<"  "<<firstfile.GetXmax()<<"  "<<firstfile.GetXmax()<<"  "<<power<<endl;
		//}
		//cout<<s1<<"  "<<power<<endl;
        //cout<<search_time<<endl;
		}}
	}
