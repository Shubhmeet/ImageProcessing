#include "opencv2/highgui.hpp"
#include <iostream>
using namespace cv;
using namespace std;

int main(int argc, char** argv) {
	if (argc != 3) {
		cout << argv[0] << ": "
			<< "got " << argc - 1 << " arguments. Expecting two: width height."
			<< endl;
		return(-1);
	}

	int width = atoi(argv[1]);
	int height = atoi(argv[2]);
	int** RED1 = new int*[height];
	int** GREEN1 = new int*[height];
	int** BLUE1 = new int*[height];
	int** RED2 = new int*[height];
	int** GREEN2 = new int*[height];
	int** BLUE2 = new int*[height];

	for (int i = 0; i < height; i++) {
		RED1[i] = new int[width];
		GREEN1[i] = new int[width];
		BLUE1[i] = new int[width];
		RED2[i] = new int[width];
		GREEN2[i] = new int[width];
		BLUE2[i] = new int[width];
	}

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
		{
			int r1, g1, b1;
			int r2, g2, b2;

			double x = (double)j / (double)width;
			double y = (double)i / (double)height;
			double Y = 1.0;

			double L = 90;
			double u = x * 512 - 255;
			double v = y * 512 - 255;


			/* Your code should be placed here
			It should translate xyY to byte sRGB
			and Luv to byte sRGB
			*/double X1, Y1, Z1;

			// checking x,y are between 0 and 1
			if (x<0) { x = 0; }
			if (x>1) { x = 1; }
			if (y<0) { y = 0; }
			if (y>1) { y = 1; }


			// conversion of xyY to XYZ-- XYZ are taken as X1Y1Z1
			X1 = (x*Y) / y;
			Y1 = Y;
			Z1 = ((1 - x - y)*Y) / y;

			double Rsrgbl = (3.24*X1) - (1.53*Y1) - (0.498*Z1);
			double Gsrgbl = (1.87*Y1) - (0.969*X1) + (0.041*Z1);
			double Bsrgbl = (0.055*X1) - (0.204*Y1) + (1.057*Z1);


			//clipping out of range values [0,1]
			if (Rsrgbl<0) {

				Rsrgbl = 0;
			}

			if (Rsrgbl>1) {
				Rsrgbl = 1;
			}



			if (Gsrgbl<0) {

				Gsrgbl = 0;
			}

			if (Gsrgbl>1) {
				Gsrgbl = 1;
			}



			if (Bsrgbl<0) {

				Bsrgbl = 0;
			}

			if (Bsrgbl>1) {
				Bsrgbl = 1;
			}


			double Rsrgbnl = 0.0;
			double Gsrgbnl = 0.0;
			double Bsrgbnl = 0.0;
			//converting srgb linear value to srgb non-linear values
			if (Rsrgbl<0.00304) {
				Rsrgbnl = 12.92*Rsrgbl;
			}
			else {
				Rsrgbnl = (1.055)*pow(Rsrgbl, 1 / 2.4) - 0.055;
			}

			if (Gsrgbl<0.00304) {
				Gsrgbnl = 12.92*Gsrgbl;
			}
			else {
				Gsrgbnl = (1.055)*pow(Gsrgbl, 1 / 2.4) - 0.055;
			}
			if (Bsrgbl<0.00304) {
				Bsrgbnl = 12.92*Bsrgbl;
			}
			else {
				Bsrgbnl = (1.055)*pow(Bsrgbl, 1 / 2.4) - 0.055;
			}

			//now streching the values of non-liner srgb values from [0,1] range to [0,255]

			r1 = (int)(Rsrgbnl * 255);
			g1 = (int)(Gsrgbnl * 255);
			b1 = (int)(Bsrgbnl * 255);



			RED1[i][j] = r1;
			GREEN1[i][j] = g1;
			BLUE1[i][j] = b1;




			double xw = 0.95;
			double yw = 1.0;
			double zw = 1.09;

			double uw = (4 * xw) / (xw + (15 * yw) + (3 * zw));

			double vw = (9 * yw) / (xw + (15 * yw) + (3 * zw));
			double u1, v1;

			if (L == 0)
			{
				u1 = 0;
				v1 = 0;

			}

			else {
				u1 = (u + (13 * uw*L)) / (13 * L);
				v1 = (v + (13 * vw*L)) / (13 * L);
			}

			double Ynew, Xnew, Znew;

			Ynew = pow(((L + 16) / 116), 3)*yw;
			if (v1 == 0)
			{
				Xnew = 0;
				Znew = 0;
			}
			else
			{
				Xnew = (Ynew*2.25*u1) / v1;
				Znew = (Ynew*(3 - (0.75*u1) - (5 * v1))) / v1;
			}

			//converting to linear SRGB values
			double linearsr, linearsg, linearsb;
			linearsr = 3.240479*Xnew - 1.53715*Ynew - 0.498535*Znew;
			linearsg = 1.875991*Ynew - 0.969256*Xnew + 0.041556*Znew;
			linearsb = 0.055648*Xnew - 0.204043*Ynew + 1.057311*Znew;




			//clippping the values linear srgb in (0-1 range)

			if (linearsr<0)
			{
				linearsr = 0;
			}

			if (linearsr>1)
			{
				linearsr = 1;
			}

			if (linearsg<0)
			{
				linearsg = 0;
			}

			if (linearsg>1)
			{
				linearsg = 1;
			}

			if (linearsb<0)
			{
				linearsb = 0;
			}

			if (linearsb>1)
			{
				linearsb = 1;
			}


			// converting linear srgb values to non linear srgb values (using gamma correction)

			double nonlinearsr, nonlinearsg, nonlinearsb;
			if (linearsr<0.00304)
			{
				nonlinearsr = 12.92*linearsr;
			}
			else
			{
				nonlinearsr = 1.055*pow(linearsr, (1 / 2.4)) - 0.055;
			}

			if (linearsg<0.00304)
			{
				nonlinearsg = 12.92*linearsg;
			}
			else
			{
				nonlinearsg = 1.055*pow(linearsg, (1 / 2.4)) - 0.055;
			}
			if (linearsb<0.00304)
			{
				nonlinearsb = 12.92*linearsb;
			}
			else
			{
				nonlinearsb = 1.055*pow(linearsb, (1 / 2.4)) - 0.055;
			}


			r2 = (int)(nonlinearsr * 255);
			g2 = (int)(nonlinearsg * 255);
			b2 = (int)(nonlinearsb * 255);


			

			// this is the end of your code

			RED1[i][j] = r1;
			GREEN1[i][j] = g1;
			BLUE1[i][j] = b1;
			RED2[i][j] = r2;
			GREEN2[i][j] = g2;
			BLUE2[i][j] = b2;
		}


	Mat R1(height, width, CV_8UC1);
	Mat G1(height, width, CV_8UC1);
	Mat B1(height, width, CV_8UC1);

	Mat R2(height, width, CV_8UC1);
	Mat G2(height, width, CV_8UC1);
	Mat B2(height, width, CV_8UC1);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++) {

			R1.at<uchar>(i, j) = RED1[i][j];
			G1.at<uchar>(i, j) = GREEN1[i][j];
			B1.at<uchar>(i, j) = BLUE1[i][j];

			R2.at<uchar>(i, j) = RED2[i][j];
			G2.at<uchar>(i, j) = GREEN2[i][j];
			B2.at<uchar>(i, j) = BLUE2[i][j];
		}

	Mat xyY;
	Mat xyY_planes[] = { B1, G1, R1 };
	merge(xyY_planes, 3, xyY);
	namedWindow("xyY", CV_WINDOW_AUTOSIZE);
	imshow("xyY", xyY);

	Mat Luv;
	Mat Luv_planes[] = { B2, G2, R2 };
	merge(Luv_planes, 3, Luv);
	namedWindow("Luv", CV_WINDOW_AUTOSIZE);
	imshow("Luv", Luv);
	imwrite("Luv.png", Luv);
	imwrite("xyY.png", xyY);
	waitKey(0); // Wait for a keystroke
	return(0);
}