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

	double** r = new double*[rows];
	double** g = new double*[rows];
	double** b = new double*[rows];

	double** x = new double*[rows];
	double** y = new double*[rows];
	double** z = new double*[rows];

	double** L = new double*[rows];
	double** u = new double*[rows];
	double** v = new double*[rows];

	double** out1 = new double*[rows];
	double** out2 = new double*[rows];
	double** out3 = new double*[rows];

	double** cr = new double*[rows];
	double** cg = new double*[rows];
	double** cb = new double*[rows];

	double** gcr = new double*[rows];
	double** gcg = new double*[rows];
	double** gcb = new double*[rows];

	int** cr1 = new int*[rows];
	int** cg1 = new int*[rows];
	int** cb1 = new int*[rows];

	double** u2 = new double*[rows];
	double** v2 = new double*[rows];

	double** lx = new double*[rows];
	double** ly = new double*[rows];
	double** lz = new double*[rows];

	double** u_dash = new double*[rows];
	double** v_dash = new double*[rows];

	int** Ld = new int*[rows];

	for (int i = 0; i < rows; i++)
	{
		R[i] = new double[cols];
		G[i] = new double[cols];
		B[i] = new double[cols];

		R1[i] = new double[cols];
		G1[i] = new double[cols];
		B1[i] = new double[cols];

		r[i] = new double[cols];
		g[i] = new double[cols];
		b[i] = new double[cols];

		x[i] = new double[cols];
		y[i] = new double[cols];
		z[i] = new double[cols];

		L[i] = new double[cols];
		u[i] = new double[cols];
		v[i] = new double[cols];

		out1[i] = new double[cols];
		out2[i] = new double[cols];
		out3[i] = new double[cols];

		cr[i] = new double[cols];
		cg[i] = new double[cols];
		cb[i] = new double[cols];

		cr1[i] = new int[cols];
		cg1[i] = new int[cols];
		cb1[i] = new int[cols];

		gcr[i] = new double[cols];
		gcg[i] = new double[cols];
		gcb[i] = new double[cols];

		u2[i] = new double[cols];
		v2[i] = new double[cols];

		lx[i] = new double[cols];
		ly[i] = new double[cols];
		lz[i] = new double[cols];

		u_dash[i] = new double[cols];
		v_dash[i] = new double[cols];

		Ld[i] = new int[cols];
	}



	//reading the input color image

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			r[i][j] = iR.at<uchar>(i, j);
			g[i][j] = iG.at<uchar>(i, j);
			b[i][j] = iB.at<uchar>(i, j);

			// convert the above image to non linear rgb

			R[i][j] = r[i][j] / 255;
			G[i][j] = g[i][j] / 255;
			B[i][j] = b[i][j] / 255;


			// Applying inverse gamma to convert non linear rgb to linear rgb-

			if (R[i][j]<0.03928)
			{
				R1[i][j] = R[i][j] / 12.92;
			}

			else
			{
				R1[i][j] = pow(((R[i][j] + 0.055) / 1.055), 2.4);
			}

			if (G[i][j]<0.03928)
			{
				G1[i][j] = G[i][j] / 12.92;
			}

			else
			{
				G1[i][j] = pow(((G[i][j] + 0.055) / 1.055), 2.4);
			}

			if (B[i][j]<0.03928)
			{
				B1[i][j] = B[i][j] / 12.92;
			}

			else
			{
				B1[i][j] = pow(((B[i][j] + 0.055) / 1.055), 2.4);
			}


			//converting to XYZ
			x[i][j] = 0.412453*R1[i][j] + 0.35758*G1[i][j] + 0.180423*B1[i][j];
			y[i][j] = 0.212671*R1[i][j] + 0.71516*G1[i][j] + 0.072169*B1[i][j];
			z[i][j] = 0.019334*R1[i][j] + 0.119193*G1[i][j] + 0.950227*B1[i][j];

			if (x[i][j]<0) x[i][j] = 0;
			if (y[i][j]<0) y[i][j] = 0;
			if (z[i][j]<0) z[i][j] = 0;


			//converting XYZ to luv
			double xw = 0.95;
			double yw = 1.0;
			double zw = 1.09;

			double uw = (4.0*xw) / (xw + 15.0*yw + 3.0*zw);
			double vw = (9.0*yw) / (xw + 15.0*yw + 3.0*zw);

			double t = y[i][j] / yw;

			if (t>0.008856)
			{
				L[i][j] = 116.0*pow(t, (1.0 / 3.0)) - 16.0;
			}

			else
			{
				L[i][j] = 903.3*t;
			}

			if (L[i][j]<0)
			{
				L[i][j] = 0.0;
			}

			if (L[i][j]>100.0)
			{
				L[i][j] = 100.0;
			}

			double d, u1, v1;
			d = x[i][j] + 15.0*y[i][j] + 3.0*z[i][j];
			u1 = (4.0*x[i][j]) / d;
			v1 = (9.0*y[i][j]) / d;
			u[i][j] = 13.0*L[i][j] * (u1 - uw);
			v[i][j] = 13.0*L[i][j] * (v1 - vw);

		}

	for (int i = 0; i<rows; i++)
		for (int j = 0; j<cols; j++)
		{
			//discretize L values
			Ld[i][j] = (int)(L[i][j] + 0.5);
		}


	//calculating the minimum and maximum L value from window
	double minl, maxl;
	minl = Ld[H1][W1];
	maxl = Ld[H1][W1];

	for (int i = H1; i <= H2; i++)
		for (int j = W1; j <= W2; j++)
		{
			if (Ld[i][j] < minl)
			{
				minl = Ld[i][j];
			}

			if (Ld[i][j] > maxl)
			{
				maxl = Ld[i][j];
			}
		}


	//Calculating the histogram
	int hist[101];
	int f[101];
	int h[101];

	for (int k = 0; k < 101; k++)
	{
		hist[k] = 0;
		f[k] = 0;
	}

	for (int i = H1; i <= H2; i++)
		for (int j = W1; j <= W2; j++)
		{

			hist[Ld[i][j]]++;
		}


	int sum = 0;
	f[0] = hist[0];
	for (int i = 1; i<101; i++)
	{
		f[i] = f[i - 1] + hist[i];

	}

	h[0] = (int)((((double)f[0] * 101.0) / (2.0* (double)f[100])));

	for (int i = 1; i<101; i++)
	{

		h[i] = (int)(((double)(f[i - 1] + f[i])*101.0) / (2.0* (double)f[100]));

	}


	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{

			if (Ld[i][j] <= minl)
			{
				Ld[i][j] = 0;
			}

			else
				if (Ld[i][j] >= maxl)
				{
					Ld[i][j] = 100;
				}

				else
					Ld[i][j] = h[Ld[i][j]];

		}


	//convert to srgb
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			double xw = 0.95;
			double yw = 1.0;
			double zw = 1.09;

			double uw = (4.0*xw) / (xw + 15.0*yw + 3.0*zw);
			double vw = (9.0*yw) / (xw + 15.0*yw + 3.0*zw);




			// converting luv to XYZ
			u_dash[i][j] = 0.0; v_dash[i][j] = 0.0;
			if (Ld[i][j]>0)
			{
				u_dash[i][j] = (u[i][j] + 13.0*uw*Ld[i][j]) / (13.0*Ld[i][j]);
				v_dash[i][j] = (v[i][j] + 13.0*vw*Ld[i][j]) / (13.0*Ld[i][j]);
			}


			if (Ld[i][j]>7.9996)
			{
				ly[i][j] = pow(((Ld[i][j] + 16.0) / 116.0), 3.0)*yw;
			}

			else
			{
				ly[i][j] = (Ld[i][j] * yw) / 903.3;
			}

			if (v_dash[i][j] == 0)
			{
				lx[i][j] = 0;
				lz[i][j] = 0;
			}

			else
			{
				lx[i][j] = (ly[i][j] * 2.25*u_dash[i][j]) / v_dash[i][j];
				lz[i][j] = (ly[i][j] * (3.0 - 0.75*u_dash[i][j] - 5.0*v_dash[i][j])) / v_dash[i][j];
			}

			//clipping XYZ values
			if (lx<0) lx = 0;
			if (ly<0) ly = 0;
			if (lz<0) lz = 0;



			//convert to SRGB

			cr[i][j] = 3.240479*lx[i][j] - 1.53715*ly[i][j] - 0.498535*lz[i][j];
			cg[i][j] = 1.875991*ly[i][j] - 0.969256*lx[i][j] + 0.041556*lz[i][j];
			cb[i][j] = 0.055648*lx[i][j] - 0.204043*ly[i][j] + 1.057311*lz[i][j];

			//clippping the values (0-1 range)

			if (cr[i][j]<0)
			{
				cr[i][j] = 0;
			}

			if (cr[i][j]>1)
			{
				cr[i][j] = 1;
			}

			if (cb[i][j]<0)
			{
				cb[i][j] = 0;
			}

			if (cb[i][j]>1)
			{
				cb[i][j] = 1;
			}

			if (cg[i][j]<0)
			{
				cg[i][j] = 0;
			}

			if (cg[i][j]>1)
			{
				cg[i][j] = 1;
			}


			// applying gamma correction to convert to non linear srgb

			if (cr[i][j]<0.00304)
			{
				gcr[i][j] = 12.92*cr[i][j];
			}

			else
			{
				gcr[i][j] = 1.055*pow(cr[i][j], (1 / 2.4)) - 0.055;
			}

			if (cg[i][j]<0.00304)
			{
				gcg[i][j] = 12.92*cg[i][j];
			}

			else
			{
				gcg[i][j] = 1.055*pow(cg[i][j], (1 / 2.4)) - 0.055;
			}

			if (cb[i][j]<0.00304)
			{
				gcb[i][j] = 12.92*cb[i][j];
			}

			else
			{
				gcb[i][j] = 1.055*pow(cb[i][j], (1 / 2.4)) - 0.055;
			}


			//stretching to 0-255 range
			cr1[i][j] = (int)((gcr[i][j] * 255) + 0.5);
			cg1[i][j] = (int)((gcg[i][j] * 255) + 0.5);
			cb1[i][j] = (int)((gcb[i][j] * 255) + 0.5);

			out1[i][j] = cr1[i][j];
			out2[i][j] = cg1[i][j];
			out3[i][j] = cb1[i][j];
		}

	Mat lr(rows, cols, CV_8UC1);
	Mat lg(rows, cols, CV_8UC1);
	Mat lb(rows, cols, CV_8UC1);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			lr.at<uchar>(i, j) = out1[i][j];;
			lg.at<uchar>(i, j) = out2[i][j];;
			lb.at<uchar>(i, j) = out3[i][j];;
		}

	Mat out_planes[] = { lb, lg, lr };
	Mat image;
	merge(out_planes, 3, image);

	namedWindow("final", CV_WINDOW_AUTOSIZE);
	imshow("final", image);
	imwrite(outName, image);

}//run

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
