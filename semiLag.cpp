/* Most of the time particle_grid_Len is a constant 100 */

#include "semiLag.hpp"

#ifdef DEBUG

#include <string>
using namespace std;

static const double interpolate(const Vec3 &v, double data[4][4][4], const Vec3 boundingBox[4][4][4]);
static const double interpolate(const Vec3 &v, double data[4][4][4], const Vec3 boundingBox[4][4][4], string s) {
	double result = interpolate(v, data, boundingBox);
	if(result < 0){
		cout << "target: " << endl;
		cout << v.dx << ' ' << v.dy << ' ' << v.dz << endl;
		cout << "boundingBox: " << endl;
		cout <<boundingBox[1][1][1].dx << ' ' << boundingBox[1][1][1].dy << ' '
			<< boundingBox[1][1][1].dz << endl;
		cout <<boundingBox[2][2][2].dx << ' ' << boundingBox[2][2][2].dy << ' '
			<< boundingBox[2][2][2].dz << endl;
	}
	return result;
}
#endif


/* p[i,j,k] at position [i * GRID_SIZE * GRID_SIZE + j * GRID_SIZE + k] */
static double pFirst[GRID_SIZE * GRID_SIZE * GRID_SIZE];
static double pSecond[GRID_SIZE * GRID_SIZE * GRID_SIZE];

static inline double dist(const Vec3 &v1, const Vec3 &v2) {
	return (v1.dx - v2.dx) * (v1.dx - v2.dx) +
		(v1.dy - v2.dy) * (v1.dy - v2.dy) +
		(v1.dz - v2.dz) * (v1.dz - v2.dz);
}

static Vec3 binary_search(const Grid &g, const Vec3 left,
		const Vec3 right, double STEP) {
	Vec3 mid = (left + right) * 0.5;
	if(dist(left, right) < STEP * STEP)
		return left;
	int midX = floor((1.0 + mid.dx) * 0.5 * GRID_SIZE);
	int midY = floor((1.0 + mid.dy) * 0.5 * GRID_SIZE);
	int midZ = floor((1.0 + mid.dz) * 0.5 * GRID_SIZE);

	if(g.particle[midX][midY][midZ].isOccupied)
		return binary_search(g, mid, right, STEP);
	else
		return binary_search(g, left, mid, STEP);
	
	/* Make gcc happy */
	return Vec3();
}

static double *poissonSolver(const Grid &g) {
	/* User of this function should provide in g 
	 * the heuristic guessing of the p value
	 * and should guarantee that v equals to u*
	 * before calculation
	 */
	const double STEP = 2.0 / GRID_SIZE;
	const double d = STEP * STEP;
	const double a = d / 6.0;
	const double af = d * d / 6.0;

	double *pold = pFirst;
	memset(pFirst, 0, GRID_SIZE * GRID_SIZE * GRID_SIZE);
	double *pnew = pSecond;
	int iter = 0;
	double maxDIFF = 0;
	do{
		/* Iterate time can be controlled for better performance */
		for(int i = 0; i != GRID_SIZE; ++i) {
			int iMinus = i - 1 < 0 ? 0 : i - 1;
			int iPlus = i + 1 >= GRID_SIZE ? i : i + 1;
			for(int j = 0; j != GRID_SIZE; ++j) {
				int jMinus = j - 1 < 0 ? 0 : j - 1;
				int jPlus = j + 1 >= GRID_SIZE ? j : j + 1;
				for(int k = 0; k != GRID_SIZE; ++k) {
					int kMinus = k - 1 < 0 ? 0 : k - 1;
					int kPlus = k + 1 >= GRID_SIZE ? k : k + 1;

					/* Here we calculate divUStar */
					double divUStar = g.particle[i][j][k].vx + 
						g.particle[i][j][k].vy + 
						g.particle[i][j][k].vz;

					/* This might be put outside the function to 
					 * care for vicinity */
					if(i - 1 >= 0)
						divUStar -= g.particle[i - 1][j][k].vx;
					if(j - 1 >= 0)
						divUStar -= g.particle[i][j - 1][k].vy;
					if(k - 1 >= 0)
						divUStar -= g.particle[i][j][k - 1].vz;

					pnew[FETCH(i, j, k)] = a * (
							pold[FETCH(iMinus, j, k)] + pold[FETCH(iPlus, j, k)] +
							pold[FETCH(i, jMinus, k)] + pold[FETCH(i, jPlus, k)] +
							pold[FETCH(i, j, kMinus)] + pold[FETCH(i, j, kPlus)]
							) - af * divUStar / TIME_STEP;
					if(pnew[FETCH(i, j, k)] - pold[FETCH(i, j, k)] > maxDIFF)
						maxDIFF = pnew[FETCH(i, j, k)] - pold[FETCH(i, j, k)];
					/*
					 *std::cout << pnew[FETCH(i,j,k)] - pold[FETCH(i,j,k)] << std::endl;
					 */
				}
			}
		}
		double *temp = pold;
		pold = pnew;
		pnew = temp;
		++iter;
	}while(iter < POISSON_ITER /*&& maxDIFF > 0.000001*/);

#ifdef DEBUG
	for(int j = 0; j != GRID_SIZE; ++j) {
		if(j != 5)
			continue;
		for(int i = 0; i != GRID_SIZE; ++i) {
			for(int k = 0; k != GRID_SIZE; ++k) {
				cout << pold[i * GRID_SIZE * GRID_SIZE + j * GRID_SIZE + k] << ' ';
			}
			cout << endl;
		}
		cout << endl;
	}
#endif
	return pold;
}


/* Hermite Interpolation */
/* This instanciation is only used for non-velocity interpolation */
static const double interpolate(const Vec3 &v, double data[4][4][4], const Vec3 boundingBox[4][4][4]) {
	/* For the first part we need to 
	 * get data from the vertex box
	 * containing the required point
	 */

	/* I think they should have been value-initialized */
	/* Actually caller of this method should 
	 * guarantee that this is true.
	 */

	/*
	 *for(int i = 0; i != 4; ++i) {
	 *    for(int j = 0; j != 4; ++j) {
	 *        for(int k = 0; k != 4; ++k) {
	 *            std::cout << boundingBox[i][j][k].dx << ' ' << boundingBox[i][j][k].dy
	 *                << ' ' << boundingBox[i][j][k].dz << std::endl;
	 *        }
	 *    }
	 *}
	 */
/*
 *    std::cout << boundingBox[1][1][1].dx << ' ' << boundingBox[1][1][1].dy
 *        << ' ' << boundingBox[1][1][1].dz << std::endl;
 *    std::cout << v.dx << ' ' << v.dy << ' ' << v.dz << ' ' << std::endl;
 *
 *    std::cout << std::endl;
 *    std::cout << std::endl;
 */

	double STEP = 2.0 / GRID_SIZE;

	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			/* Interpolate along third axis */
			double diff = (v.dz - boundingBox[i][j][1].dz) / STEP;

			if(std::abs(diff) < 0.0001) {
				data[i][j][0] = data[i][j][1];
				
				/*
				 *if(data[i][j][0] < 0)
				 *    std::cout << boundingBox[i][j][1].dz << std::endl;
				 */

				continue;
			}
			if(std::abs(diff - 1.0) < 0.0001) {
				data[i][j][0] = data[i][j][2];
				continue;
			}
			double quadDiff = diff * diff;
			double cubicDiff = quadDiff * diff;
			double deltaK = data[i][j][2] - data[i][j][1];
			double dk = (data[i][j][2] - data[i][j][0]) * 0.5;
			double dkPlus = (data[i][j][3] - data[i][j][1]) * 0.5;

			/*
			 *dk = (dk * deltaK < 0 ? 0 : dk);
			 *dkPlus = (dkPlus * deltaK < 0 ? 0 : dkPlus);
			 */
			/*
			 *if(dk * deltaK <= 0)
			 *    dk = 0;
			 *if(dkPlus * deltaK <= 0)
			 *    dkPlus = 0;
			 */

			/* The article got this wrong */
			double coeffA3 = dk + dkPlus - deltaK * 2;
			double coeffA2 = deltaK * 3 - dk * 2 - dkPlus;
			double coeffA1 = dk;
			double coeffA0 = data[i][j][1];

			if(coeffA3 * deltaK <= 0 || coeffA2 * deltaK <= 0 || coeffA1 * deltaK <= 0) {
				coeffA3 = - deltaK * 2;
				coeffA2 = deltaK * 3;
				coeffA1 = 0;
			}

			data[i][j][0] = coeffA3 * cubicDiff + coeffA2 * quadDiff + coeffA1 * diff + coeffA0;
		}
		/* Interpolate along second axis */
		double diff = (v.dy - boundingBox[i][1][1].dy) / STEP;
		if(std::abs(diff) < 0.0001) {
			data[i][0][0] = data[i][1][0];

			/*
			 *if(data[i][0][0] < 0)
			 *    std::cout << data[i][0][0] << std::endl;
			 */

			continue;
		}
		if(std::abs(diff - 1.0) < 0.0001) {
			data[i][0][0] = data[i][2][0];
			continue;
		}
		double quadDiff = diff * diff;
		double cubicDiff = quadDiff * diff;
		double deltaK = data[i][2][0] - data[i][1][0];
		double dk = (data[i][2][0] - data[i][0][0]) * 0.5;
		double dkPlus = (data[i][3][0] - data[i][1][0]) * 0.5;

		/*
		 *dk = (dk * deltaK < 0 ? 0 : dk);
		 *dkPlus = (dkPlus * deltaK < 0 ? 0 : dkPlus);
		 */
		/*
		 *if(dk * deltaK <= 0)
		 *    dk = 0;
		 *if(dkPlus * deltaK <= 0)
		 *    dkPlus = 0;
		 */

		/* The article got this wrong */
		double coeffA3 = dk + dkPlus - deltaK * 2;
		double coeffA2 = deltaK * 3 - dk * 2 - dkPlus;
		double coeffA1 = dk;
		double coeffA0 = data[i][1][0];

		if(coeffA3 * deltaK <= 0 || coeffA2 * deltaK <= 0 || coeffA1 * deltaK <= 0) {
			coeffA3 = - deltaK * 2;
			coeffA2 = deltaK * 3;
			coeffA1 = 0;
		}

		data[i][0][0] = coeffA3 * cubicDiff + coeffA2 * quadDiff + 
			coeffA1 * diff + coeffA0;
		/*
		 *if(std::abs(data[i][0][0] + 0.0308606) < 0.001) {
		 *    std::cout <<"dk:" << dk << " dkPlus" << dkPlus << " deltaK" << deltaK
		 *        << std::endl;
		 *    std::cout <<coeffA3 << ' ' << coeffA2 << ' ' << coeffA1 << ' '
		 *        << coeffA0 << std::endl;
		 *    std::cout << diff << std::endl;
		 *}
		 */
	}
	/* Interpolate along third axis */
	double diff = (v.dx - boundingBox[1][1][1].dx) / STEP;
	if(std::abs(diff) < 0.0001) {

		/*
		 *if(data[0][0][0] < 0)
		 *    std::cout << data[0][0][0] << std::endl;
		 */

		return data[1][0][0];
	}
	if(std::abs(diff - 1.0) < 0.0001) {
		return data[2][0][0];
	}
	double quadDiff = diff * diff;
	double cubicDiff = quadDiff * diff;
	double deltaK = data[2][0][0] - data[1][0][0];
	double dk = (data[2][0][0] - data[0][0][0]) * 0.5;
	double dkPlus = (data[3][0][0] - data[1][0][0]) * 0.5;

	/*
	 *dk = (dk * deltaK < 0 ? 0 : dk);
	 *dkPlus = (dkPlus * deltaK < 0 ? 0 : dkPlus);
	 */
	/*
	 *if(dk * deltaK <= 0)
	 *    dk = 0;
	 *if(dkPlus * deltaK <= 0)
	 *    dkPlus = 0;
	 */

	/* The article got this wrong */
	double coeffA3 = dk + dkPlus - deltaK * 2;
	double coeffA2 = deltaK * 3 - dk * 2 - dkPlus;
	double coeffA1 = dk;
	double coeffA0 = data[1][0][0];

	if(coeffA3 * deltaK <= 0 || coeffA2 * deltaK <= 0 || coeffA1 * deltaK <= 0) {
		coeffA3 = - deltaK * 2;
		coeffA2 = deltaK * 3;
		coeffA1 = 0;
	}

	data[0][0][0] = coeffA3 * cubicDiff + coeffA2 * quadDiff + coeffA1 * diff + coeffA0;

	/*
	 *if(data[0][0][0] < 0)
	 *    std::cout << data[0][0][0] << std::endl;
	 */

	return data[0][0][0];
}

/* Now we have to instantiate different interpolationWrapper */
static Vec2 tempXdensYInterpolationWrapper(const Vec3&loc, const Grid &g) {
	int xLow = floor((1.0 + loc.dx) * 0.5 * GRID_SIZE - 0.5);
	int yLow = floor((1.0 + loc.dy) * 0.5 * GRID_SIZE - 0.5);
	int zLow = floor((1.0 + loc.dz) * 0.5 * GRID_SIZE - 0.5);
	/*
	 *std::cout << loc.dx << ' ' << loc.dy << ' ' << loc.dz << std::endl;
	 *std::cout << xLow << ' ' << yLow << ' ' << zLow << std::endl;
	 */
	double STEP = 2.0 / GRID_SIZE;
	double dataTemp[4][4][4];
	double dataDens[4][4][4];
	Vec3 bBox[4][4][4];
	memset(dataTemp, 0.0, sizeof(dataTemp));
	memset(dataDens, 0.0, sizeof(dataDens));
	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			for(int k = 0; k != 4; ++k) {
				bBox[i][j][k].dx = (xLow + i - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dy = (yLow + j - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dz = (zLow + k - 1 + 0.5) * STEP - 1.0;
				if(xLow + i - 1 < 0 || xLow + i - 1 >= GRID_SIZE)
					continue;
				if(yLow + j - 1 < 0 || yLow + j - 1 >= GRID_SIZE)
					continue;
				if(zLow + k - 1 < 0 || zLow + k - 1 >= GRID_SIZE)
					continue;
				dataTemp[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].temperature;
				dataDens[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].density;
				/*
				 *if(dataDens[i][j][k] < 0 )
				 *    std::cout << "less" << std::endl;
				 */
			}
		}
	}

	/*
	 *if(bBox[1][1][1].dz - loc.dz > 0) {
	 *    std::cout << zLow << ' ' << loc.dz << std::endl;
	 *}
	 */

	Vec2 v;
	v.dx = interpolate(loc, dataTemp, bBox);
	v.dy = interpolate(loc, dataDens, bBox);

	return v;
}

static Vec3 velocityInterpolationWrapper(const Vec3 &loc, const Grid &g) {
	Vec3 result;
	int xLow = floor((1.0 + loc.dx) * 0.5 * GRID_SIZE - 1.0);
	int yLow = floor((1.0 + loc.dy) * 0.5 * GRID_SIZE - 0.5);
	int zLow = floor((1.0 + loc.dz) * 0.5 * GRID_SIZE - 0.5);
	double STEP = 2.0 / GRID_SIZE;
	double dataVX[4][4][4];
	Vec3 bBox[4][4][4];
	memset(dataVX, 0, sizeof(dataVX));
	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			for(int k = 0; k != 4; ++k) {
				bBox[i][j][k].dx = (xLow + i - 1 + 1.0) * STEP - 1.0;
				bBox[i][j][k].dy = (yLow + j - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dz = (zLow + k - 1 + 0.5) * STEP - 1.0;
				if(xLow + i - 1 < 0 || xLow + i - 1 >= GRID_SIZE)
					continue;
				if(yLow + j - 1 < 0 || yLow + j - 1 >= GRID_SIZE)
					continue;
				if(zLow + k - 1 < 0 || zLow + k - 1 >= GRID_SIZE)
					continue;
				dataVX[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vStarX;
			}
		}
	}
	result.dx = interpolate(loc, dataVX, bBox);
	xLow = floor((1.0 + loc.dx) * 0.5 * GRID_SIZE - 0.5);
	yLow = floor((1.0 + loc.dy) * 0.5 * GRID_SIZE - 1.0);
	zLow = floor((1.0 + loc.dz) * 0.5 * GRID_SIZE - 0.5);
	double dataVY[4][4][4];
	memset(dataVY, 0, sizeof(dataVY));
	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			for(int k = 0; k != 4; ++k) {
				bBox[i][j][k].dx = (xLow + i - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dy = (yLow + j - 1 + 1.0) * STEP - 1.0;
				bBox[i][j][k].dz = (zLow + k - 1 + 0.5) * STEP - 1.0;
				if(xLow + i - 1 < 0 || xLow + i - 1 >= GRID_SIZE)
					continue;
				if(yLow + j - 1 < 0 || yLow + j - 1 >= GRID_SIZE)
					continue;
				if(zLow + k - 1 < 0 || zLow + k - 1 >= GRID_SIZE)
					continue;
				dataVY[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vStarY;
			}
		}
	}

	/*
	 *if(bBox[1][1][1].dz - loc.dz > 0) {
	 *    std::cout << zLow << ' ' << loc.dz << std::endl;
	 *}
	 */

	result.dy = interpolate(loc, dataVY, bBox);
	xLow = floor((1.0 + loc.dx) * 0.5 * GRID_SIZE - 0.5);
	yLow = floor((1.0 + loc.dy) * 0.5 * GRID_SIZE - 0.5);
	zLow = floor((1.0 + loc.dz) * 0.5 * GRID_SIZE - 1.0);
	double dataVZ[4][4][4];
	memset(dataVZ, 0, sizeof(dataVZ));
	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			for(int k = 0; k != 4; ++k) {
				bBox[i][j][k].dx = (xLow + i - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dy = (yLow + j - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dz = (zLow + k - 1 + 1.0) * STEP - 1.0;
				if(xLow + i - 1 < 0 || xLow + i - 1 >= GRID_SIZE)
					continue;
				if(yLow + j - 1 < 0 || yLow + j - 1 >= GRID_SIZE)
					continue;
				if(zLow + k - 1 < 0 || zLow + k - 1 >= GRID_SIZE)
					continue;
				dataVZ[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vStarZ;
			}
		}
	}
	result.dz = interpolate(loc, dataVZ, bBox);

	/*
	 *std::cout << result.dz << std::endl;
	 */

	return result;
}

static double vXInterpolateWrapper(const Vec3 &loc, const Grid &g) {
	int xLow = floor((1.0 + loc.dx) * 0.5 * GRID_SIZE - 1.0);
	int yLow = floor((1.0 + loc.dy) * 0.5 * GRID_SIZE - 0.5);
	int zLow = floor((1.0 + loc.dz) * 0.5 * GRID_SIZE - 0.5);
	double STEP = 2.0 / GRID_SIZE;
	double dataVX[4][4][4];
	Vec3 bBox[4][4][4];
	memset(dataVX, 0, sizeof(dataVX));
	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			for(int k = 0; k != 4; ++k) {
				bBox[i][j][k].dx = (xLow + i - 1 + 1.0) * STEP - 1.0;
				bBox[i][j][k].dy = (yLow + j - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dz = (zLow + k - 1 + 0.5) * STEP - 1.0;
				if(xLow + i - 1 < 0 || xLow + i - 1 >= GRID_SIZE)
					continue;
				if(yLow + j - 1 < 0 || yLow + j - 1 >= GRID_SIZE)
					continue;
				if(zLow + k - 1 < 0 || zLow + k - 1 >= GRID_SIZE)
					continue;
				dataVX[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vx;
			}
		}
	}
	return interpolate(loc, dataVX, bBox);
}

static double vYInterpolateWrapper(const Vec3 &loc, const Grid &g) {
	int xLow = floor((1.0 + loc.dx) * 0.5 * GRID_SIZE - 0.5);
	int yLow = floor((1.0 + loc.dy) * 0.5 * GRID_SIZE - 1.0);
	int zLow = floor((1.0 + loc.dz) * 0.5 * GRID_SIZE - 0.5);
	double STEP = 2.0 / GRID_SIZE;
	double dataVY[4][4][4];
	Vec3 bBox[4][4][4];
	memset(dataVY, 0, sizeof(dataVY));
	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			for(int k = 0; k != 4; ++k) {
				bBox[i][j][k].dx = (xLow + i - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dy = (yLow + j - 1 + 1.0) * STEP - 1.0;
				bBox[i][j][k].dz = (zLow + k - 1 + 0.5) * STEP - 1.0;
				if(xLow + i - 1 < 0 || xLow + i - 1 >= GRID_SIZE)
					continue;
				if(yLow + j - 1 < 0 || yLow + j - 1 >= GRID_SIZE)
					continue;
				if(zLow + k - 1 < 0 || zLow + k - 1 >= GRID_SIZE)
					continue;
				dataVY[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vy;
			}
		}
	}
	return interpolate(loc, dataVY, bBox);
}

static double vZInterpolateWrapper(const Vec3 &loc, const Grid &g) {
	int xLow = floor((1.0 + loc.dx) * 0.5 * GRID_SIZE - 0.5);
	int yLow = floor((1.0 + loc.dy) * 0.5 * GRID_SIZE - 0.5);
	int zLow = floor((1.0 + loc.dz) * 0.5 * GRID_SIZE - 1.0);
	double STEP = 2.0 / GRID_SIZE;
	double dataVZ[4][4][4];
	Vec3 bBox[4][4][4];
	memset(dataVZ, 0, sizeof(dataVZ));
	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 4; ++j) {
			for(int k = 0; k != 4; ++k) {
				bBox[i][j][k].dx = (xLow + i - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dy = (yLow + j - 1 + 0.5) * STEP - 1.0;
				bBox[i][j][k].dz = (zLow + k - 1 + 1.0) * STEP - 1.0;
				if(xLow + i - 1 < 0 || xLow + i - 1 >= GRID_SIZE)
					continue;
				if(yLow + j - 1 < 0 || yLow + j - 1 >= GRID_SIZE)
					continue;
				if(zLow + k - 1 < 0 || zLow + k - 1 >= GRID_SIZE)
					continue;
				dataVZ[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vz;
			}
		}
	}
	return interpolate(loc, dataVZ, bBox);
}

/* 
 * We need to calculate advection in multiple dimension(3d)
 * and for temperature and density the equation looks like
 * dF/dt = ∂F/∂t + [u(x,t) * ∆] ∂F/∂x
 * for this we need the result F(x - 2a, t - ∆t)
 */

/* This calculator presumes that the velocity field of the 
 * last step i.e. that in the old Grid is one with forcing 
 * terms added back */

/* UPDATE 2017/07/02 */
void semiLagrangeCalc(const Grid &old, Grid &gen) {
	/* Then we'll solve the advection term */
	double STEP = 2.0 / GRID_SIZE;


	for(int i = 0; i != GRID_SIZE; ++i) {
		for(int j = 0; j != GRID_SIZE; ++j) {
			for(int k = 0; k != GRID_SIZE; ++k) {
				Vec3 curr(
						-1 + STEP * (i + 0.5),
						-1 + STEP * (j + 0.5),
						-1 + STEP * (k + 0.5)
						);
				/* Solve equation use original advectionTerm */
				Vec3 alpha = Vec3(old.particle[i][j][k].vStarX,
						old.particle[i][j][k].vStarY,
						old.particle[i][j][k].vStarZ) * TIME_STEP * 0.5;
				Vec3 newAlpha = velocityInterpolationWrapper(curr - alpha * 0.5, old)
					* TIME_STEP;
				int iter = 1;
				while(dist(alpha, newAlpha) >= CONTROL && iter < ITER_TIMES) {
					alpha = newAlpha;
					/* use vStar to interpolate */
					newAlpha = velocityInterpolationWrapper(curr - alpha * 0.5, old)
						* TIME_STEP;
					++iter;
				}
				gen.particle[i][j][k].advectionTerm = newAlpha;
			}
		}
	}


	/* Then according to that term, T, D, vNew
	 * will all be interpolated */

	for(int i = 0; i != GRID_SIZE; ++i) {
		for(int j = 0; j != GRID_SIZE; ++j) {
			for(int k = 0; k != GRID_SIZE; ++k) {
				
				/* This part can be varied if new movement option is activated */
				gen.particle[i][j][k].boarder = old.particle[i][j][k].boarder;
				gen.particle[i][j][k].isOccupied = old.particle[i][j][k].isOccupied;
				/* We can use determined moving velocity to calculate new object position */

				if(gen.particle[i][j][k].isOccupied) {
					gen.particle[i][j][k].vx = old.particle[i][j][k].vx;
					gen.particle[i][j][k].vy = old.particle[i][j][k].vy;
					gen.particle[i][j][k].vz = old.particle[i][j][k].vz;
					continue;
				}

				Vec3 curr(
						-1 + STEP * (i + 0.5),
						-1 + STEP * (j + 0.5),
						-1 + STEP * (k + 0.5)
						);
				Vec3 target = curr - gen.particle[i][j][k].advectionTerm;

				int midX = floor((1.0 + target.dx) * 0.5 * GRID_SIZE);
				int midY = floor((1.0 + target.dy) * 0.5 * GRID_SIZE);
				int midZ = floor((1.0 + target.dz) * 0.5 * GRID_SIZE);

				if(
						midX >= 0 && midX < GRID_SIZE &&
						midY >= 0 && midY < GRID_SIZE &&
						midZ >= 0 && midZ < GRID_SIZE &&
						old.particle[midX][midY][midZ].isOccupied
						) {
					target = binary_search(old, target, curr, 0.10 * STEP);
					midX = floor((1.0 + target.dx) * 0.5 * GRID_SIZE);
					midY = floor((1.0 + target.dy) * 0.5 * GRID_SIZE);
					midZ = floor((1.0 + target.dz) * 0.5 * GRID_SIZE);

					gen.particle[i][j][k].vx = old.particle[midX][midY][midZ].vx;
					gen.particle[i][j][k].vy = old.particle[midX][midY][midZ].vy;
					gen.particle[i][j][k].vz = old.particle[midX][midY][midZ].vz;

					continue;
				}
				
				/* Then interpolate vX, vY, vZ */
				gen.particle[i][j][k].vx = vXInterpolateWrapper(target, old);
				gen.particle[i][j][k].vy = vYInterpolateWrapper(target, old);
				gen.particle[i][j][k].vz = vZInterpolateWrapper(target, old);
			}
		}
	}


	/* Noticed that velocity at point [i, j, k] is of center point of voxel */
	for(int i = 0; i != GRID_SIZE; ++i) {
		int iPlus = i + 1 >= GRID_SIZE ? i : i + 1;
		for(int j = 0; j != GRID_SIZE; ++j) {
			int jPlus = j + 1 >= GRID_SIZE ? j : j + 1;
			for(int k = 0; k != GRID_SIZE; ++k) {
				int kPlus = k + 1 >= GRID_SIZE ? k : k + 1;

				if(gen.particle[i][j][k].isOccupied) {
					gen.particle[i][j][k].vx = old.particle[i][j][k].vx;
					gen.particle[i][j][k].vy = old.particle[i][j][k].vy;
					gen.particle[i][j][k].vz = old.particle[i][j][k].vz;
				}
				gen.particle[i][j][k].vx = 
					(gen.particle[i][j][k].vx + gen.particle[iPlus][j][k].vx) / 2.0;
				gen.particle[i][j][k].vy = 
					(gen.particle[i][j][k].vy + gen.particle[i][jPlus][k].vy) / 2.0;
				gen.particle[i][j][k].vz = 
					(gen.particle[i][j][k].vx + gen.particle[i][j][kPlus].vx) / 2.0;
			}
		}
	}

	double *p = poissonSolver(gen);

	for(int i = 0; i != GRID_SIZE; ++i) {
		int iPlus = i + 1 >= GRID_SIZE ? i : i + 1;
		for(int j = 0; j != GRID_SIZE; ++j) {
			int jPlus = j + 1 >= GRID_SIZE ? j : j + 1;
			for(int k = 0; k != GRID_SIZE; ++k) {
				int kPlus = k + 1 >= GRID_SIZE ? k : k + 1;
				/* First subtract influence of p to 
				 * make the fluid incompressible
				 */
				if(gen.particle[i][j][k].isOccupied)
					continue;
				double dpX = iPlus != i ?
					(p[FETCH(iPlus, j, k)] - p[FETCH(i, j, k)]) / STEP : 0;
				double dpY = jPlus != j ?
					(p[FETCH(i, jPlus, k)] - p[FETCH(i, j, k)]) / STEP : 0;
				double dpZ = kPlus != k ?
					(p[FETCH(i, j, kPlus)] - p[FETCH(i, j, k)]) / STEP : 0;
				gen.particle[i][j][k].vx -= dpX * TIME_STEP;
				gen.particle[i][j][k].vy -= dpY * TIME_STEP;
				gen.particle[i][j][k].vz -= dpZ * TIME_STEP;
#ifdef DEBUG
				/*
				 *cout << dpX << ' ' << dpY << ' ' << dpZ << endl;
				 */
#endif
				/* We'll calculate vStar next */
				gen.particle[i][j][k].vStarX = 1.5 * gen.particle[i][j][k].vx
					- 0.5 * old.particle[i][j][k].vx;
				gen.particle[i][j][k].vStarY = 1.5 * gen.particle[i][j][k].vy
					- 0.5 * old.particle[i][j][k].vy;
				gen.particle[i][j][k].vStarZ = 1.5 * old.particle[i][j][k].vz
					- 0.5 * old.particle[i][j][k].vz;

				/*
				 *std::cout << old.particle[i][j][k].vx
				 *    << ' ' << old.particle[i][j][k].vy
				 *    << ' ' << old.particle[i][j][k].vz << std::endl;
				 */

			}
		}
	}

	for(int i = 0; i != GRID_SIZE; ++i) {
		for(int j = 0; j != GRID_SIZE; ++j) {
			for(int k = 0; k != GRID_SIZE; ++k) {
				Vec3 curr(
						-1 + STEP * (i + 0.5),
						-1 + STEP * (j + 0.5),
						-1 + STEP * (k + 0.5)
						);
				/* Solve equation use original advectionTerm */
				Vec3 alpha = Vec3(gen.particle[i][j][k].vStarX,
						old.particle[i][j][k].vStarY,
						old.particle[i][j][k].vStarZ) * TIME_STEP * 0.5;
				Vec3 newAlpha = velocityInterpolationWrapper(curr - alpha * 0.5, gen)
					* TIME_STEP;
				int iter = 1;
				while(dist(alpha, newAlpha) >= CONTROL && iter < ITER_TIMES) {
					alpha = newAlpha;
					/* use vStar to interpolate */
					newAlpha = velocityInterpolationWrapper(curr - alpha * 0.5, gen)
						* TIME_STEP;
					++iter;
				}
				gen.particle[i][j][k].advectionTerm = newAlpha;
			}
		}
	}

	for(int i = 0; i != GRID_SIZE; ++i) {
		for(int j = 0; j != GRID_SIZE; ++j) {
			for(int k = 0; k != GRID_SIZE; ++k) {
				if(gen.particle[i][j][k].isOccupied)
					continue;

				Vec3 curr(
						-1 + STEP * (i + 0.5),
						-1 + STEP * (j + 0.5),
						-1 + STEP * (k + 0.5)
						);
				Vec3 target = curr - gen.particle[i][j][k].advectionTerm;

				int midX = floor((1.0 + target.dx) * 0.5 * GRID_SIZE);
				int midY = floor((1.0 + target.dy) * 0.5 * GRID_SIZE);
				int midZ = floor((1.0 + target.dz) * 0.5 * GRID_SIZE);

				Vec2 tempXdensY;

				if(
						midX >= 0 && midX < GRID_SIZE &&
						midY >= 0 && midY < GRID_SIZE &&
						midZ >= 0 && midZ < GRID_SIZE &&
						old.particle[midX][midY][midZ].isOccupied
						) {
					target = binary_search(old, target, curr, 0.10 * STEP);
					tempXdensY = tempXdensYInterpolationWrapper(target, old);

					gen.particle[i][j][k].temperature = tempXdensY.dx;
					gen.particle[i][j][k].density = tempXdensY.dy;

					gen.particle[i][j][k].advectionTerm = target;

					continue;
				}
				tempXdensY = tempXdensYInterpolationWrapper(target, old);

				gen.particle[i][j][k].temperature = tempXdensY.dx;
				gen.particle[i][j][k].density = tempXdensY.dy;

			}
		}
	}

	for(int i = 0; i != GRID_SIZE; ++i) {
		for(int j = 0; j != GRID_SIZE; ++j) {
			for(int k = 0; k != GRID_SIZE; ++k) {
				if(gen.particle[i][j][k].isOccupied) {
					/* Grid that is occupied is made zero_density */
					if(gen.particle[i][j][k].boarder) {
						int count = 0;
						if(i - 1 >= 0 && !gen.particle[i - 1][j][k].isOccupied) {
							gen.particle[i][j][k].density =
								gen.particle[i - 1][j][k].density;
							++count;
						}
						if(j - 1 >= 0 && !gen.particle[i][j - 1][k].isOccupied) {
							gen.particle[i][j][k].density +=
								gen.particle[i][j - 1][k].density;
							++count;
						}
						if(k - 1 >= 0 && !gen.particle[i][j][k - 1].isOccupied) {
							gen.particle[i][j][k].density +=
								gen.particle[i][j][k - 1].density;
							++count;
						}
						if(i + 1 < GRID_SIZE && !gen.particle[i + 1][j][k].isOccupied) {
							gen.particle[i][j][k].density +=
								gen.particle[i + 1][j][k].density;
							++count;
						}
						if(j + 1 < GRID_SIZE && !gen.particle[i][j + 1][k].isOccupied) {
							gen.particle[i][j][k].density +=
								gen.particle[i][j + 1][k].density;
							++count;
						}
						if(k + 1 < GRID_SIZE && !gen.particle[i][j][k + 1].isOccupied) {
							gen.particle[i][j][k].density +=
								gen.particle[i][j][k + 1].density;
							++count;
						}
						gen.particle[i][j][k].density /= static_cast<double>(count);
					}
					else gen.particle[i][j][k].density = 0.0;
				}
			}
		}
	}

}
