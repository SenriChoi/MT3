#include <Novice.h>
#include <cmath>

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


Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio,float nearClip, float farClip ) {
	Matrix4x4 result = {};
	result.m[0][0] = (1 / aspectRatio) * (1 / std::tan(fovY / 2));
	result.m[1][1] = (1 / std::tan(fovY / 2));
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1;
	result.m[3][2] = (- 1 * nearClip* farClip)/ (farClip - nearClip);

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

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height,float minDepth, float maxDepth) {
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

static const int kRowHeight = 20;
static const int kColumnWidth = 60;


void MatrixScreenPrintf(int x, int y, const Matrix4x4& matrix, const char* name) {
	Novice::ScreenPrintf(x,y,"%s",name);
	for (int row = 1; row <= 4; ++row) {
		for (int column = 0; column <4; ++column) {
				Novice::ScreenPrintf(x + column * kColumnWidth,
					y + row * kRowHeight, "%6.02f", matrix.m[row - 1][column ]);
			
		}
	}
};


// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {


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
		Matrix4x4 orthographicMatrix = MakeOrthographicMatrix(-160.0f,160.0f,200.0f,300.0f,0.0f,1000.0f);
		Matrix4x4 perspectiveFovMatrix = MakePerspectiveFovMatrix(0.63f,1.33f,0.1f,1000.0f);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(100.0f,200.0f,600.0f,300.0f,0.0f,1.0f);

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		MatrixScreenPrintf(0, 0, orthographicMatrix, "orthographicMatrix");
		MatrixScreenPrintf(0, kRowHeight*5, perspectiveFovMatrix, "perspectiveFovMatrix");
		MatrixScreenPrintf(0, kRowHeight * 5*2, viewportMatrix, "viewportMatrix");


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
