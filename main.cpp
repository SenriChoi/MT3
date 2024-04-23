#include <Novice.h>
#include <cmath>
#include <assert.h>
const char kWindowTitle[] = "GC2A_06_ジョ_カエイ";

struct Matrix4x4 {
	float m[4][4];
};

struct Vector3 {
	float x;
	float y;
	float z;
};

//積
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

#pragma region 逆行列
//逆行列
Matrix4x4 Inverse(const Matrix4x4& matrix) {
	Matrix4x4 result{};
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

Vector3 Cross(const Vector3& v1, const Vector3& v2) {
	Vector3 result = {};
	result.x = v1.y * v2.z - v1.z * v2.y;
	result.y = v1.z * v2.x - v1.x * v2.z;
	result.z = v1.x * v2.y - v1.y* v2.x;

	return result;
};

static const int kRowHeight = 20;
static const int kColumnWidth = 60;


void VectorScreenPrintf(int x, int y, const Vector3& vector3, const char* name) {
	Novice::ScreenPrintf(x, y, "%6.02f  %6.02f  %6.02f  %s", vector3.x, vector3.y, vector3.z, name);
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {
	Vector3 v1{ 1.2f,-3.9f,2.5f };
	Vector3 v2{ 2.8f,0.4f,-1.3f };
	Vector3 rotate{ 0,0,0 };
	Vector3 translate{ 0,0,0 };
	Vector3 cameraPosition{ 0,0,-20 };
	Vector3 kLocalVertices[3] = { {1,-1,0},{0,1,0},{-1,-1,0} };

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///
		/// 
	
		rotate.y += 0.05f;
		if (rotate.y >= 3.14f*2) {

			rotate.y = 0;
		}

		if (keys[DIK_W]) {
			translate.z -= 0.1f;

		}
		else if (keys[DIK_S]) {
			translate.z += 0.1f;
		}


		if (keys[DIK_A]) {
			translate.x -= 0.1f;

		}
		else if (keys[DIK_D]) {
			translate.x += 0.1f;
		}

		Vector3 cross = Cross(v1, v2);

	

		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, cameraPosition);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, (1280.0f/720.0f), 0.1f, 100.0f);
		Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0.0f,0.0f,1280.0f,720.0f,0.0f,1.0f);
		Vector3 screenVertices[3];
		for (uint32_t i = 0; i < 3; ++i) {
			Vector3 ndcVertex = Transform(kLocalVertices[i], worldViewProjectionMatrix);
			screenVertices[i] = Transform(ndcVertex, viewportMatrix);
		}


		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		VectorScreenPrintf(0, 0, cross, "Cross");
		Novice::DrawTriangle((int)screenVertices[0].x, (int)screenVertices[0].y, (int)screenVertices[1].x, (int)screenVertices[1].y, (int)screenVertices[2].x, (int)screenVertices[2].y,RED,kFillModeSolid);
		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
