#include "opencv2/highgui.hpp"
#include <iostream>
using namespace cv;
using namespace std;

void runOnWindow(int W1, int H1, int W2, int H2, Mat inputImage, char *outName) {
	int rows = inputImage.rows;
	int cols = inputImage.cols;

	vector<Mat> i_planes;
	split(inputImage, i_planes);
	Mat iB = i_planes[0];
	Mat iG = i_planes[1];
	Mat iR = i_planes[2];

	// dynamically allocate RGB arrays of size rows x cols
	double** R = new double*[rows];
	double** G = new double*[rows];
	double** B = new double*[rows];

	double** R1 = new double*[rows];
	double** G1 = new double*[rows];
	double** B1 = new double*[rows];


	double** Ri = new double*[rows];
	double** Gi = new double*[rows];
	double** Bi = new double*[rows];

	double** X = new double*[rows];
	double** Y = new double*[rows];
	double** Z = new double*[rows];

	double** x = new double*[rows];
	double** y = new double*[rows];

	double** Xnew = new double*[rows];
	double** Ynew = new double*[rows];
	double** Znew = new double*[rows];


	double** RED1 = new double*[rows];
	double** GREEN1 = new double*[rows];
	double** BLUE1 = new double*[rows];

	double** linearsr = new double*[rows];
	double** linearsg = new double*[rows];
	double** linearsb = new double*[rows];

	double** nonlinearsr = new double*[rows];
	double** nonlinearsg = new double*[rows];
	double** nonlinearsb = new double*[rows];



	for (int i = 0; i < rows; i++) {
		R[i] = new double[cols];
		G[i] = new double[cols];
		B[i] = new double[cols];
		R1[i] = new double[cols];
		G1[i] = new double[cols];
		B1[i] = new double[cols];
		Ri[i] = new double[cols];
		Gi[i] = new double[cols];
		Bi[i] = new double[cols];

		X[i] = new double[cols];
		Y[i] = new double[cols];
		Z[i] = new double[cols];

		x[i] = new double[cols];
		y[i] = new double[cols];

		Xnew[i] = new double[cols];

		Ynew[i] = new double[cols];
		Znew[i] = new double[cols];



		RED1[i] = new double[cols];
		GREEN1[i] = new double[cols];
		BLUE1[i] = new double[cols];

		linearsr[i] = new double[cols];
		linearsg[i] = new double[cols];
		linearsb[i] = new double[cols];


		nonlinearsr[i] = new double[cols];
		nonlinearsg[i] = new double[cols];
		nonlinearsb[i] = new double[cols];


	}

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {
			R[i][j] = iR.at<uchar>(i, j);
			G[i][j] = iG.at<uchar>(i, j);
			B[i][j] = iB.at<uchar>(i, j);


			//converting nonlinear srgb in[0,255] range to non linear[0,1]



			R1[i][j] = R[i][j] / 255;
			G1[i][j] = G[i][j] / 255;
			B1[i][j] = B[i][j] / 255;

			if (R1[i][j]<0.03928)
			{
				Ri[i][j] = R1[i][j] / 12.92;
			}
			else
			{
				Ri[i][j] = pow((R1[i][j] + 0.055) / 1.055, 2.4);
			}


			if (G1[i][j]<0.03928)
			{
				Gi[i][j] = G1[i][j] / 12.92;
			}
			else
			{
				Gi[i][j] = pow((G1[i][j] + 0.055) / 1.055, 2.4);
			}


			if (B1[i][j]<0.03928)
			{
				Bi[i][j] = B1[i][j] / 12.92;
			}
			else
			{
				Bi[i][j] = pow((B1[i][j] + 0.055) / 1.055, 2.4);
			}



			X[i][j] = (0.412453*Ri[i][j]) + (0.35758*Gi[i][j]) + (0.180423*Bi[i][j]);
			Y[i][j] = (0.212671*Ri[i][j]) + (0.71516*Gi[i][j]) + (0.072169*Bi[i][j]);
			Z[i][j] = (0.019334*Ri[i][j]) + (0.119193*Gi[i][j]) + (0.950227*Bi[i][j]);

			if ((X[i][j] + Y[i][j] + Z[i][j]) == 0) {
				x[i][j] = 0;
				y[i][j] = 0;

			}
			else {

				x[i][j] = X[i][j] / (X[i][j] + Y[i][j] + Z[i][j]);
				y[i][j] = Y[i][j] / (X[i][j] + Y[i][j] + Z[i][j]);
			}


		}


	//	   The transformation should be based on the
	//	   historgram of the pixels in the W1,W2,H1,H2 range.
	//	   The following code goes over these pixels
	double minL = 200;
	double maxL = 0;
	for (int i = H1; i <= H2; i++)
		for (int j = W1; j <= W2; j++) {
			if (Y[i][j]<minL) {

				minL = Y[i][j];

			}

			if (Y[i][j]>maxL) {

				maxL = Y[i][j];

			}


			/* double r = R[i][j];
			double g = G[i][j];
			double b = B[i][j];
			int gray = (int) (0.3*r + 0.6*g + 0.1*b + 0.5);

			R[i][j] = G[i][j] = B[i][j] = gray;*/
		}


	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {
			/*	if(Y[i][j]<=minL){

			Y[i][j]=0;

			}

			if(Y[i][j]>=maxL){

			Y[i][j]=1;

			}*/



			if ((maxL - minL) == 0) {
				Y[i][j] = 0;
			}
			else {
				Y[i][j] = (Y[i][j] - minL)*(1) / (maxL - minL);
			}




		}


	for (int i = 0; i<rows; i++)
		for (int j = 0; j<cols; j++) {


			if (y[i][j] == 0) {

				Xnew[i][j] = 0;
				Znew[i][j] = 0;

			}
			else
			{
				Xnew[i][j] = (x[i][j] * Y[i][j]) / y[i][j];
				Znew[i][j] = ((1 - x[i][j] - y[i][j])*Y[i][j]) / y[i][j];
				Ynew[i][j] = Y[i][j];
			}


			//converting to linear SRGB values

			linearsr[i][j] = 3.240479*Xnew[i][j] - 1.53715*Ynew[i][j] - 0.498535*Znew[i][j];
			linearsg[i][j] = 1.875991*Ynew[i][j] - 0.969256*Xnew[i][j] + 0.041556*Znew[i][j];
			linearsb[i][j] = 0.055648*Xnew[i][j] - 0.204043*Ynew[i][j] + 1.057311*Znew[i][j];




			//clippping the values linear srgb in (0-1 range)

			if (linearsr[i][j]<0)
			{
				linearsr[i][j] = 0;
			}

			if (linearsr[i][j]>1)
			{
				linearsr[i][j] = 1;
			}

			if (linearsg[i][j]<0)
			{
				linearsg[i][j] = 0;
			}

			if (linearsg[i][j]>1)
			{
				linearsg[i][j] = 1;
			}

			if (linearsb[i][j]<0)
			{
				linearsb[i][j] = 0;
			}

			if (linearsb[i][j]>1)
			{
				linearsb[i][j] = 1;
			}


			// converting linear srgb values to non linear srgb values (using gamma correction)

			if (linearsr[i][j]<0.00304)
			{
				nonlinearsr[i][j] = 12.92*linearsr[i][j];
			}
			else
			{
				nonlinearsr[i][j] = 1.055*pow(linearsr[i][j], (1 / 2.4)) - 0.055;
			}

			if (linearsg[i][j]<0.00304)
			{
				nonlinearsg[i][j] = 12.92*linearsg[i][j];
			}
			else
			{
				nonlinearsg[i][j] = 1.055*pow(linearsg[i][j], (1 / 2.4)) - 0.055;
			}
			if (linearsb[i][j]<0.00304)
			{
				nonlinearsb[i][j] = 12.92*linearsb[i][j];
			}
			else
			{
				nonlinearsb[i][j] = 1.055*pow(linearsb[i][j], (1 / 2.4)) - 0.055;
			}

			int finalsr = (int)(nonlinearsr[i][j] * 255);
			int finalsg = (int)(nonlinearsg[i][j] * 255);
			int finalsb = (int)(nonlinearsb[i][j] * 255);

			RED1[i][j] = finalsr;
			GREEN1[i][j] = finalsg;
			BLUE1[i][j] = finalsb;


		}


	Mat oR(rows, cols, CV_8UC1);
	Mat oG(rows, cols, CV_8UC1);
	Mat oB(rows, cols, CV_8UC1);


	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {


			oR.at<uchar>(i, j) = RED1[i][j];;
			oG.at<uchar>(i, j) = GREEN1[i][j];;
			oB.at<uchar>(i, j) = BLUE1[i][j];;


		}
	Mat o_planes[] = { oB, oG, oR };
	Mat outImage;
	merge(o_planes, 3, outImage);

	namedWindow("output", CV_WINDOW_AUTOSIZE);
	//
	imshow("output", outImage);
	imwrite(outName, outImage);
}

int main(int argc, char** argv) {
	if (argc != 7) {
		cerr << argv[0] << ": "
			<< "got " << argc - 1
			<< " arguments. Expecting six: w1 h1 w2 h2 ImageIn ImageOut."
			<< endl;
		cerr << "Example: proj1b 0.2 0.1 0.8 0.5 fruits.jpg out.bmp" << endl;
		return(-1);
	}
	double w1 = atof(argv[1]);
	double h1 = atof(argv[2]);
	double w2 = atof(argv[3]);
	double h2 = atof(argv[4]);
	char *inputName = argv[5];
	char *outputName = argv[6];

	if (w1<0 || h1<0 || w2 <= w1 || h2 <= h1 || w2>1 || h2>1) {
		cerr << " arguments must satisfy 0 <= w1 < w2 <= 1"
			<< " ,  0 <= h1 < h2 <= 1" << endl;
		return(-1);
	}

	Mat inputImage = imread(inputName, CV_LOAD_IMAGE_UNCHANGED);
	if (inputImage.empty()) {
		cout << "Could not open or find the image " << inputName << endl;
		return(-1);
	}

	string windowInput("input: ");
	windowInput += inputName;

	namedWindow(windowInput, CV_WINDOW_AUTOSIZE);
	imshow(windowInput, inputImage);

	if (inputImage.type() != CV_8UC3) {
		cout << inputName << " is not a standard color image  " << endl;
		return(-1);
	}

	int rows = inputImage.rows;
	int cols = inputImage.cols;
	int W1 = (int)(w1*(cols - 1));
	int H1 = (int)(h1*(rows - 1));
	int W2 = (int)(w2*(cols - 1));
	int H2 = (int)(h2*(rows - 1));

	runOnWindow(W1, H1, W2, H2, inputImage, outputName);

	waitKey(0); // Wait for a keystroke
	return(0);
}
