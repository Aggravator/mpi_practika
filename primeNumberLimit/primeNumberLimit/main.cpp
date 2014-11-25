#include <stdio.h>
#include <time.h>
#include <conio.h>
const int LIMIT=100000000;
int main(){
	int i;
	bool step=false;
	clock_t start,end;
	start=clock();
	for(i=0;i<LIMIT;i+=step?1:3){
		step=!step;
	}
	end=clock();
	printf("step over=%ld\n",end-start);
	start=clock();
	for(i=0;i<LIMIT;++i);
	end=clock();
	printf("line=%ld\n",end-start);
	start=clock();
	for(i=0;i<LIMIT;i+=4);
	end=clock();
	printf("2 gap=%ld\n",end-start);
	_getch();
}