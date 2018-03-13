#include <stdio.h>

//rotate/flip a quadrant appropriately
void rot(int s, int *x, int *y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            printf("%d, %d\n",*x, *y);
            *x = n-1 - *x;
            *y = n-1 - *y;
            printf("%d, %d\n",*x, *y);
        }

        //Swap x and y
        int t  = *x;
        *x = *y;
        *y = t;
    }
}

//convert d to (x,y)
void d2xy(int n, int d, int *x, int *y) {
    int rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s*=2) {
        rx = 1 & (t/2);
        ry = 1 & (t ^ rx);
        printf("%d, %d, %d, %d, %d, %d\n", s, *x, *y, rx, ry, t);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}


int main(){

	int x, y;
	int n = 8;
	int d = 3;
	d2xy(n, d, &x, &y);
	printf("The position is (%d, %d)\n", x, y);
	return 0;
}
