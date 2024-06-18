#pragma once
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

struct Vector3 {
	float x;
	float y;
	float z;
};

struct Matrix4x4 {
	float m[4][4];
};

Vector3 Subtract(const Vector3& v1, const Vector3& v2) {

	return{ v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };

};
Vector3 Multiply(float scaler, const Vector3& v) {
	Vector3 result = { 0,0,0 };
	result.x = scaler * v.x;
	result.y = scaler * v.y;
	result.z = scaler * v.z;
	return result;
}
Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result = { 0,0,0 };
	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;
	return result;
}
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result;

	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];


	assert(w != 0.0f);

	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
};
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return result;
};
#pragma region MakeAffineMatrix

Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
	Matrix4x4 result = {};

	result.m[0][0] = 1.0f;
	result.m[1][1] = 1.0f;
	result.m[2][2] = 1.0f;
	result.m[3][3] = 1.0f;
	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;

	return result;
};

Matrix4x4 MakeScaleMatrix(const Vector3& scale) {
	Matrix4x4 result = {};

	result.m[0][0] = scale.x;
	result.m[1][1] = scale.y;
	result.m[2][2] = scale.z;
	result.m[3][3] = 1.0f;

	return result;
};

Matrix4x4 MakeRotateXMatrix(float radian) {
	Matrix4x4 result = {};
	result.m[0][0] = 1;
	result.m[3][3] = 1;
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = std::sin(radian);
	result.m[2][1] = -1 * std::sin(radian);
	result.m[2][2] = std::cos(radian);

	return result;
};
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result = {};
	result.m[1][1] = 1;
	result.m[3][3] = 1;
	result.m[0][0] = std::cos(radian);
	result.m[0][2] = -1 * std::sin(radian);
	result.m[2][0] = std::sin(radian);
	result.m[2][2] = std::cos(radian);

	return result;

};
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result = {};
	result.m[2][2] = 1;
	result.m[3][3] = 1;
	result.m[0][0] = std::cos(radian);
	result.m[0][1] = std::sin(radian);
	result.m[1][0] = -1 * std::sin(radian);
	result.m[1][1] = std::cos(radian);

	return result;

};


Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 result = {};
	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);

	Matrix4x4 rotateXYZMatrix = Multiply(rotateXMatrix, Multiply(rotateYMatrix, rotateZMatrix));

	result = Multiply(Multiply(MakeScaleMatrix(scale), rotateXYZMatrix), MakeTranslateMatrix(translate));


	return result;
};

#pragma endregion
#pragma region 逆行列
//逆行列
Matrix4x4 Inverse(const Matrix4x4& matrix) {
	Matrix4x4 result = {};
	float determinant = matrix.m[0][0] * (matrix.m[1][1] * matrix.m[2][2] * matrix.m[3][3] +
		matrix.m[2][1] * matrix.m[3][2] * matrix.m[1][3] +
		matrix.m[3][1] * matrix.m[1][2] * matrix.m[2][3] -
		matrix.m[3][1] * matrix.m[2][2] * matrix.m[1][3] -
		matrix.m[2][1] * matrix.m[1][2] * matrix.m[3][3] -
		matrix.m[1][1] * matrix.m[3][2] * matrix.m[2][3]) -
		matrix.m[0][1] * (matrix.m[1][0] * matrix.m[2][2] * matrix.m[3][3] +
			matrix.m[2][0] * matrix.m[3][2] * matrix.m[1][3] +
			matrix.m[3][0] * matrix.m[1][2] * matrix.m[2][3] -
			matrix.m[3][0] * matrix.m[2][2] * matrix.m[1][3] -
			matrix.m[2][0] * matrix.m[1][2] * matrix.m[3][3] -
			matrix.m[1][0] * matrix.m[3][2] * matrix.m[2][3]) +
		matrix.m[0][2] * (matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][3] +
			matrix.m[2][0] * matrix.m[3][1] * matrix.m[1][3] +
			matrix.m[3][0] * matrix.m[1][1] * matrix.m[2][3] -
			matrix.m[3][0] * matrix.m[2][1] * matrix.m[1][3] -
			matrix.m[2][0] * matrix.m[1][1] * matrix.m[3][3] -
			matrix.m[1][0] * matrix.m[3][1] * matrix.m[2][3]) -
		matrix.m[0][3] * (matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][2] +
			matrix.m[2][0] * matrix.m[3][1] * matrix.m[1][2] +
			matrix.m[3][0] * matrix.m[1][1] * matrix.m[2][2] -
			matrix.m[3][0] * matrix.m[2][1] * matrix.m[1][2] -
			matrix.m[2][0] * matrix.m[1][1] * matrix.m[3][2] -
			matrix.m[1][0] * matrix.m[3][1] * matrix.m[2][2]);



	if (determinant != 0) {
		result.m[0][0] = (matrix.m[1][1] * matrix.m[2][2] * matrix.m[3][3] +
			matrix.m[2][1] * matrix.m[3][2] * matrix.m[1][3] +
			matrix.m[3][1] * matrix.m[1][2] * matrix.m[2][3] -
			matrix.m[3][1] * matrix.m[2][2] * matrix.m[1][3] -
			matrix.m[2][1] * matrix.m[1][2] * matrix.m[3][3] -
			matrix.m[1][1] * matrix.m[3][2] * matrix.m[2][3]) /
			determinant;

		result.m[0][1] = -(matrix.m[0][1] * matrix.m[2][2] * matrix.m[3][3] +
			matrix.m[2][1] * matrix.m[3][2] * matrix.m[0][3] +
			matrix.m[3][1] * matrix.m[0][2] * matrix.m[2][3] -
			matrix.m[3][1] * matrix.m[2][2] * matrix.m[0][3] -
			matrix.m[2][1] * matrix.m[0][2] * matrix.m[3][3] -
			matrix.m[0][1] * matrix.m[3][2] * matrix.m[2][3]) /
			determinant;

		result.m[0][2] = (matrix.m[0][1] * matrix.m[1][2] * matrix.m[3][3] +
			matrix.m[1][1] * matrix.m[3][2] * matrix.m[0][3] +
			matrix.m[3][1] * matrix.m[0][2] * matrix.m[1][3] -
			matrix.m[3][1] * matrix.m[1][2] * matrix.m[0][3] -
			matrix.m[1][1] * matrix.m[0][2] * matrix.m[3][3] -
			matrix.m[0][1] * matrix.m[3][2] * matrix.m[1][3]) /
			determinant;

		result.m[0][3] = -(matrix.m[0][1] * matrix.m[1][2] * matrix.m[2][3] +
			matrix.m[1][1] * matrix.m[2][2] * matrix.m[0][3] +
			matrix.m[2][1] * matrix.m[0][2] * matrix.m[1][3] -
			matrix.m[2][1] * matrix.m[1][2] * matrix.m[0][3] -
			matrix.m[1][1] * matrix.m[0][2] * matrix.m[2][3] -
			matrix.m[0][1] * matrix.m[2][2] * matrix.m[1][3]) /
			determinant;


		result.m[1][0] = -(matrix.m[1][0] * matrix.m[2][2] * matrix.m[3][3] +
			matrix.m[2][0] * matrix.m[3][2] * matrix.m[1][3] +
			matrix.m[3][0] * matrix.m[1][2] * matrix.m[2][3] -
			matrix.m[3][0] * matrix.m[2][2] * matrix.m[1][3] -
			matrix.m[2][0] * matrix.m[1][2] * matrix.m[3][3] -
			matrix.m[1][0] * matrix.m[3][2] * matrix.m[2][3]) /
			determinant;

		result.m[1][1] = (matrix.m[0][0] * matrix.m[2][2] * matrix.m[3][3] +
			matrix.m[2][0] * matrix.m[3][2] * matrix.m[0][3] +
			matrix.m[3][0] * matrix.m[0][2] * matrix.m[2][3] -
			matrix.m[3][0] * matrix.m[2][2] * matrix.m[0][3] -
			matrix.m[2][0] * matrix.m[0][2] * matrix.m[3][3] -
			matrix.m[0][0] * matrix.m[3][2] * matrix.m[2][3]) /
			determinant;

		result.m[1][2] = -(matrix.m[0][0] * matrix.m[1][2] * matrix.m[3][3] +
			matrix.m[1][0] * matrix.m[3][2] * matrix.m[0][3] +
			matrix.m[3][0] * matrix.m[0][2] * matrix.m[1][3] -
			matrix.m[3][0] * matrix.m[1][2] * matrix.m[0][3] -
			matrix.m[1][0] * matrix.m[0][2] * matrix.m[3][3] -
			matrix.m[0][0] * matrix.m[3][2] * matrix.m[1][3]) /
			determinant;

		result.m[1][3] = (matrix.m[0][0] * matrix.m[1][2] * matrix.m[2][3] +
			matrix.m[1][0] * matrix.m[2][2] * matrix.m[0][3] +
			matrix.m[2][0] * matrix.m[0][2] * matrix.m[1][3] -
			matrix.m[2][0] * matrix.m[1][2] * matrix.m[0][3] -
			matrix.m[1][0] * matrix.m[0][2] * matrix.m[2][3] -
			matrix.m[0][0] * matrix.m[2][2] * matrix.m[1][3]) /
			determinant;


		result.m[2][0] = (matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][3] +
			matrix.m[2][0] * matrix.m[3][1] * matrix.m[1][3] +
			matrix.m[3][0] * matrix.m[1][1] * matrix.m[2][3] -
			matrix.m[3][0] * matrix.m[2][1] * matrix.m[1][3] -
			matrix.m[2][0] * matrix.m[1][1] * matrix.m[3][3] -
			matrix.m[1][0] * matrix.m[3][1] * matrix.m[2][3]) /
			determinant;

		result.m[2][1] = -(matrix.m[0][0] * matrix.m[2][1] * matrix.m[3][3] +
			matrix.m[2][0] * matrix.m[3][1] * matrix.m[0][3] +
			matrix.m[3][0] * matrix.m[0][1] * matrix.m[2][3] -
			matrix.m[3][0] * matrix.m[2][1] * matrix.m[0][3] -
			matrix.m[2][0] * matrix.m[0][1] * matrix.m[3][3] -
			matrix.m[0][0] * matrix.m[3][1] * matrix.m[2][3]) /
			determinant;

		result.m[2][2] = (matrix.m[0][0] * matrix.m[1][1] * matrix.m[3][3] +
			matrix.m[1][0] * matrix.m[3][1] * matrix.m[0][3] +
			matrix.m[3][0] * matrix.m[0][1] * matrix.m[1][3] -
			matrix.m[3][0] * matrix.m[1][1] * matrix.m[0][3] -
			matrix.m[1][0] * matrix.m[0][1] * matrix.m[3][3] -
			matrix.m[0][0] * matrix.m[3][1] * matrix.m[1][3]) /
			determinant;

		result.m[2][3] = -(matrix.m[0][0] * matrix.m[1][1] * matrix.m[2][3] +
			matrix.m[1][0] * matrix.m[2][1] * matrix.m[0][3] +
			matrix.m[2][0] * matrix.m[0][1] * matrix.m[1][3] -
			matrix.m[2][0] * matrix.m[1][1] * matrix.m[0][3] -
			matrix.m[1][0] * matrix.m[0][1] * matrix.m[2][3] -
			matrix.m[0][0] * matrix.m[2][1] * matrix.m[1][3]) /
			determinant;

		result.m[3][0] = -(matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][2] +
			matrix.m[2][0] * matrix.m[3][1] * matrix.m[1][2] +
			matrix.m[3][0] * matrix.m[1][1] * matrix.m[2][2] -
			matrix.m[3][0] * matrix.m[2][1] * matrix.m[1][2] -
			matrix.m[2][0] * matrix.m[1][1] * matrix.m[3][2] -
			matrix.m[1][0] * matrix.m[3][1] * matrix.m[2][2]) /
			determinant;

		result.m[3][1] = (matrix.m[0][0] * matrix.m[2][1] * matrix.m[3][2] +
			matrix.m[2][0] * matrix.m[3][1] * matrix.m[0][2] +
			matrix.m[3][0] * matrix.m[0][1] * matrix.m[2][2] -
			matrix.m[3][0] * matrix.m[2][1] * matrix.m[0][2] -
			matrix.m[2][0] * matrix.m[0][1] * matrix.m[3][2] -
			matrix.m[0][0] * matrix.m[3][1] * matrix.m[2][2]) /
			determinant;

		result.m[3][2] = -(matrix.m[0][0] * matrix.m[1][1] * matrix.m[3][2] +
			matrix.m[1][0] * matrix.m[3][1] * matrix.m[0][2] +
			matrix.m[3][0] * matrix.m[0][1] * matrix.m[1][2] -
			matrix.m[3][0] * matrix.m[1][1] * matrix.m[0][2] -
			matrix.m[1][0] * matrix.m[0][1] * matrix.m[3][2] -
			matrix.m[0][0] * matrix.m[3][1] * matrix.m[1][2]) /
			determinant;

		result.m[3][3] = (matrix.m[0][0] * matrix.m[1][1] * matrix.m[2][2] +
			matrix.m[1][0] * matrix.m[2][1] * matrix.m[0][2] +
			matrix.m[2][0] * matrix.m[0][1] * matrix.m[1][2] -
			matrix.m[2][0] * matrix.m[1][1] * matrix.m[0][2] -
			matrix.m[1][0] * matrix.m[0][1] * matrix.m[2][2] -
			matrix.m[0][0] * matrix.m[2][1] * matrix.m[1][2]) /
			determinant;
	}

	return result;
};
#pragma endregion
float Dot(const Vector3& v1, const Vector3& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }

float Length(const Vector3& v) {
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

float LengthSquared(const Vector3& v) {
	return Dot(v, v);
}

struct Line {
	Vector3 origin;//始
	Vector3 diff;//終
};

struct Ray {
	Vector3 origin;//始
	Vector3 diff;//終
};

struct Segment {
	Vector3 origin;//始
	Vector3 diff;//終
};


#pragma region モデル画面

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {};
	result.m[0][0] = (1 / aspectRatio) * (1 / std::tan(fovY / 2));
	result.m[1][1] = (1 / std::tan(fovY / 2));
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1;
	result.m[3][2] = (-1 * nearClip * farClip) / (farClip - nearClip);

	return result;
};

Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip) {
	Matrix4x4 result = {};
	result.m[0][0] = 2 / (right - left);
	result.m[1][1] = 2 / (top - bottom);
	result.m[2][2] = 1 / (farClip - nearClip);
	result.m[3][0] = (left + right) / (left - right);
	result.m[3][1] = (top + bottom) / (bottom - top);
	result.m[3][2] = nearClip / (nearClip - farClip);
	result.m[3][3] = 1;
	return result;
};

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 result = {};

	result.m[0][0] = width / 2;
	result.m[1][1] = -1 * (height / 2);
	result.m[2][2] = maxDepth - minDepth;
	result.m[3][0] = left + (width / 2);
	result.m[3][1] = top + (height / 2);
	result.m[3][2] = minDepth;
	result.m[3][3] = 1;

	return result;
}
#pragma endregion




//網 グリッド
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfwidth = 2.0f;
	const uint32_t kSubdivision = 10;
	const float kGridEvery = (kGridHalfwidth * 2.0f) / float(kSubdivision);

	Matrix4x4 worldViewProjectionMatrix = viewProjectionMatrix;
	Matrix4x4 viewPortMatrix = viewportMatrix;

	Vector3 startPointPosHorizontal[kSubdivision + 1];
	Vector3 endPointPosHorizontal[kSubdivision + 1];
	Vector3 startPointScreenHorizontal[kSubdivision + 1];
	Vector3 endPointScreenHorizontal[kSubdivision + 1];

	Vector3 startPointPosVertical[kSubdivision + 1];
	Vector3 endPointPosVertical[kSubdivision + 1];
	Vector3 startPointScreenVertical[kSubdivision + 1];
	Vector3 endPointScreenVertical[kSubdivision + 1];


	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		startPointPosHorizontal[xIndex].x = -kGridHalfwidth;
		endPointPosHorizontal[xIndex].x = kGridHalfwidth;

		startPointPosHorizontal[xIndex].y = 0;
		endPointPosHorizontal[xIndex].y = 0;

		startPointPosHorizontal[xIndex].z = kGridHalfwidth - kGridEvery * xIndex;
		endPointPosHorizontal[xIndex].z = kGridHalfwidth - kGridEvery * xIndex;

		Vector3 ndcStartPointPosHorizontal = Transform(startPointPosHorizontal[xIndex], worldViewProjectionMatrix);
		Vector3 ndcEndPointPosHorizontal = Transform(endPointPosHorizontal[xIndex], worldViewProjectionMatrix);

		startPointScreenHorizontal[xIndex] = Transform(ndcStartPointPosHorizontal, viewPortMatrix);;
		endPointScreenHorizontal[xIndex] = Transform(ndcEndPointPosHorizontal, viewPortMatrix);

		if (xIndex == 5) {
			Novice::DrawLine((int)startPointScreenHorizontal[xIndex].x, (int)startPointScreenHorizontal[xIndex].y, (int)endPointScreenHorizontal[xIndex].x, (int)endPointScreenHorizontal[xIndex].y, BLACK);
		}
		else {
			Novice::DrawLine((int)startPointScreenHorizontal[xIndex].x, (int)startPointScreenHorizontal[xIndex].y, (int)endPointScreenHorizontal[xIndex].x, (int)endPointScreenHorizontal[xIndex].y, WHITE);
		}

	}


	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
		startPointPosVertical[zIndex].x = kGridHalfwidth - kGridEvery * zIndex;
		endPointPosVertical[zIndex].x = kGridHalfwidth - kGridEvery * zIndex;

		startPointPosVertical[zIndex].y = 0;
		endPointPosVertical[zIndex].y = 0;

		startPointPosVertical[zIndex].z = -kGridHalfwidth;
		endPointPosVertical[zIndex].z = kGridHalfwidth;

		Vector3 ndcStartPointPosVertical = Transform(startPointPosVertical[zIndex], worldViewProjectionMatrix);
		Vector3 ndcEndPointPosVertical = Transform(endPointPosVertical[zIndex], worldViewProjectionMatrix);

		startPointScreenVertical[zIndex] = Transform(ndcStartPointPosVertical, viewPortMatrix);;
		endPointScreenVertical[zIndex] = Transform(ndcEndPointPosVertical, viewPortMatrix);

		if (zIndex == 5) {
			Novice::DrawLine((int)startPointScreenVertical[zIndex].x, (int)startPointScreenVertical[zIndex].y, (int)endPointScreenVertical[zIndex].x, (int)endPointScreenVertical[zIndex].y, BLACK);
		}
		else {
			Novice::DrawLine((int)startPointScreenVertical[zIndex].x, (int)startPointScreenVertical[zIndex].y, (int)endPointScreenVertical[zIndex].x, (int)endPointScreenVertical[zIndex].y, WHITE);
		}

	}
}


//球
struct Sphere
{
	Vector3 center;//中心点
	float radius;
};
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	//Matrix4x4 worldViewProjectionMatrix = viewProjectionMatrix;
	//Matrix4x4 viewPortMatrix = viewportMatrix;
	const uint32_t kSubdivision = 20;
	const float kLonEvery = 2.0f * float(M_PI) / float(kSubdivision);
	const float kLatEvery = float(M_PI) / float(kSubdivision);

	//Vector3 a[kSubdivision + 1];
	//Vector3 b[kSubdivision + 1];
	//Vector3 c[kSubdivision + 1];

	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -float(M_PI) / 2.0f + kLatEvery * float(latIndex);		//現在緯度

		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = kLonEvery * float(lonIndex);		//現在の経度
			Vector3 a;
			Vector3 b;
			Vector3 c;

			a.x = sphere.center.x + sphere.radius * cos(lat) * cos(lon);
			a.y = sphere.center.y + sphere.radius * sin(lat);
			a.z = sphere.center.z + sphere.radius * cos(lat) * sin(lon);

			b.x = sphere.center.x + sphere.radius * cos(lat + kLatEvery) * cos(lon);
			b.y = sphere.center.y + sphere.radius * sin(lat + kLatEvery);
			b.z = sphere.center.z + sphere.radius * cos(lat + kLatEvery) * sin(lon);

			c.x = sphere.center.x + sphere.radius * cos(lat) * cos(lon + kLonEvery);
			c.y = sphere.center.y + sphere.radius * sin(lat);
			c.z = sphere.center.z + sphere.radius * cos(lat) * sin(lon + kLonEvery);
			//座標
			Vector3 ndcA = Transform(a, viewProjectionMatrix);
			Vector3 ndcB = Transform(b, viewProjectionMatrix);
			Vector3 ndcC = Transform(c, viewProjectionMatrix);
			Vector3 screenA = Transform(ndcA, viewportMatrix);
			Vector3 screenB = Transform(ndcB, viewportMatrix);
			Vector3 screenC = Transform(ndcC, viewportMatrix);

			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x, (int)screenB.y, color);
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenC.x, (int)screenC.y, color);
		}

	}

}


//正射影ベクトルと最近接点
Vector3 Project(const Vector3& v1, const Vector3& v2) {
	float dot = Dot(v1, v2);
	float len = LengthSquared(v2);
	float t = dot / len;
	return Multiply(t, v2);
}
Vector3 Closest(const Vector3& point, const Segment& segment) {

	Vector3 pointO = segment.origin;
	Vector3 pointP = point;
	Vector3 pointA = Subtract(pointP, pointO);
	Vector3 pointB = segment.diff;


	Vector3 closestPoint = Project(pointA, pointB);

	return Add(pointO, closestPoint);
}