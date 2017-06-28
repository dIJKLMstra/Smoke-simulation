/* Most of the time particle_grid_Len is a constant 100 */

#include "semiLag.hpp"


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
			double quadDiff = diff * diff;
			double cubicDiff = quadDiff * diff;
			double deltaK = data[i][j][2] - data[i][j][1];
			double dk = (data[i][j][2] - data[i][j][0]) * 0.5;
			double dkPlus = (data[i][j][3] - data[i][j][1]) * 0.5;

			/*
			 *dk = (dk * deltaK < 0 ? 0 : dk);
			 *dkPlus = (dkPlus * deltaK < 0 ? 0 : dkPlus);
			 */
			if(dk * deltaK <= 0 || dkPlus * deltaK <= 0)
				dk = dkPlus = 0;

			/* The article got this wrong */
			double coeffA3 = dk + dkPlus - deltaK * 2;
			double coeffA2 = deltaK * 3 - dk * 2 - dkPlus;
			double coeffA1 = dk;
			double coeffA0 = data[i][j][1];
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
		double quadDiff = diff * diff;
		double cubicDiff = quadDiff * diff;
		double deltaK = data[i][2][0] - data[i][1][0];
		double dk = (data[i][2][0] - data[i][0][0]) * 0.5;
		double dkPlus = (data[i][3][0] - data[i][1][0]) * 0.5;

		/*
		 *dk = (dk * deltaK < 0 ? 0 : dk);
		 *dkPlus = (dkPlus * deltaK < 0 ? 0 : dkPlus);
		 */
		if(dk * deltaK <= 0 || dkPlus * deltaK <= 0)
			dk = dkPlus = 0;

		/* The article got this wrong */
		double coeffA3 = dk + dkPlus - deltaK * 2;
		double coeffA2 = deltaK * 3 - dk * 2 - dkPlus;
		double coeffA1 = dk;
		double coeffA0 = data[i][1][0];
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
	double quadDiff = diff * diff;
	double cubicDiff = quadDiff * diff;
	double deltaK = data[2][0][0] - data[1][0][0];
	double dk = (data[2][0][0] - data[0][0][0]) * 0.5;
	double dkPlus = (data[3][0][0] - data[1][0][0]) * 0.5;

	/*
	 *dk = (dk * deltaK < 0 ? 0 : dk);
	 *dkPlus = (dkPlus * deltaK < 0 ? 0 : dkPlus);
	 */
	if(dk * deltaK <= 0 || dkPlus * deltaK <= 0)
		dk = dkPlus = 0;

	/* The article got this wrong */
	double coeffA3 = dk + dkPlus - deltaK * 2;
	double coeffA2 = deltaK * 3 - dk * 2 - dkPlus;
	double coeffA1 = dk;
	double coeffA0 = data[1][0][0];
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

	/*
	 *if(v.dy < 0) {
	 *    std::cout << loc.dx << ' ' << loc.dy << ' ' << loc.dz << std::endl;
	 *    std::cout << bBox[1][1][1].dx << ' ' << 
	 *        bBox[1][1][1].dy << ' ' << bBox[1][1][1].dz << std::endl;
	 *}
	 */
	/*
	 *std::cout << v.dx << ' ' << v.dy << std::endl;
	 */
	/*
	 *if(v.dy < 0) {
	 *    std::cout << v.dy << std::endl;
	 *    for(int i = 0; i != 4; ++i) {
	 *        for(int j = 0; j != 4; ++j) {
	 *            for(int k = 0; k != 4; ++k) {
	 *                std::cout << dataDens[i][j][k] << ' ';
	 *            }
	 *            std::cout << std::endl;
	 *        }
	 *        std::cout << std::endl << std::endl;
	 *    }
	 *}
	 */

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
				dataVX[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vx;
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
				dataVY[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vy;
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
				dataVZ[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vz;
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
				dataVX[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vStarX;
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
				dataVY[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vStarY;
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
				dataVZ[i][j][k] = g.particle[xLow + i - 1][yLow + j - 1][zLow + k - 1].vStarZ;
			}
		}
	}
	return interpolate(loc, dataVZ, bBox);
}

static inline double dist(const Vec3 &v1, const Vec3 &v2) {
	return (v1.dx - v2.dx) * (v1.dx - v2.dx) +
		(v1.dy - v2.dy) * (v1.dy - v2.dy) +
		(v1.dz - v2.dz) * (v1.dz - v2.dz);
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
				Vec3 alpha = old.particle[i][j][k].advectionTerm;
				Vec3 newAlpha = velocityInterpolationWrapper(curr - alpha * 0.5, old)
					* TIME_STEP;
				int iter = 1;
				while(dist(alpha, newAlpha) < 0.10 && iter < ITER_TIMES) {
					alpha = newAlpha;
					/* use vStar to interpolate */
					newAlpha = velocityInterpolationWrapper(curr - alpha * 0.5, old)
						* TIME_STEP;
					++iter;
				}

				/*
				 *alpha = newAlpha;
				 *std::cout << i << ' ' << j<< ' ' << k << std::endl;
				 *std::cout << "alpha " << alpha.dx << ' ' << alpha.dy
				 *    << ' ' << alpha.dz << std::endl;
				 */

				gen.particle[i][j][k].advectionTerm = newAlpha;
			}
		}
	}
	/* Then according to that term, T, D, vNew
	 * will all be interpolated */
	for(int i = 0; i != GRID_SIZE; ++i) {
		for(int j = 0; j != GRID_SIZE; ++j) {
			for(int k = 0; k != GRID_SIZE; ++k) {
				Vec3 curr(
						-1 + STEP * (i + 0.5),
						-1 + STEP * (j + 0.5),
						-1 + STEP * (k + 0.5)
						);
				Vec2 tempXdensY = tempXdensYInterpolationWrapper(
						curr - gen.particle[i][j][k].advectionTerm, old);
				gen.particle[i][j][k].temperature = tempXdensY.dx;
				gen.particle[i][j][k].density = tempXdensY.dy;
				
				/* Then interpolate vX, vY, vZ */
				Vec3 vXAlpha = gen.particle[i][j][k].advectionTerm;
				if(i + 1 < GRID_SIZE) {
					vXAlpha = 
						(vXAlpha + gen.particle[i + 1][j][k].advectionTerm) * 0.5;
				}
				Vec3 vYAlpha = gen.particle[i][j][k].advectionTerm;
				if(j + 1 < GRID_SIZE) {
					vYAlpha = 
						(vYAlpha + gen.particle[i][j + 1][k].advectionTerm) * 0.5;
				}
				Vec3 vZAlpha = gen.particle[i][j][k].advectionTerm;
				if(j + 1 < GRID_SIZE) {
					vZAlpha = 
						(vZAlpha + gen.particle[i][j][k + 1].advectionTerm) * 0.5;
				}
				Vec3 locX = curr;
				locX.dx += 0.5 * STEP;
				Vec3 locY = curr;
				locY.dy += 0.5 * STEP;
				Vec3 locZ = curr;
				locZ.dz += 0.5 * STEP;
				gen.particle[i][j][k].vx = vXInterpolateWrapper(locX - vXAlpha, old);
				gen.particle[i][j][k].vy = vYInterpolateWrapper(locY - vYAlpha, old);
				gen.particle[i][j][k].vz = vZInterpolateWrapper(locZ - vZAlpha, old);

				/* Then we will update V* */
				gen.particle[i][j][k].vStarX = 1.5 * gen.particle[i][j][k].vx
					- 0.5 * old.particle[i][j][k].vx;
				gen.particle[i][j][k].vStarY = 1.5 * gen.particle[i][j][k].vy
					- 0.5 * old.particle[i][j][k].vy;
				gen.particle[i][j][k].vStarZ = 1.5 * gen.particle[i][j][k].vz
					- 0.5 * old.particle[i][j][k].vz;
			}
		}
	}
}
