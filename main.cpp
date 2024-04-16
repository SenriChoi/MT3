#include <Novice.h>
#include <math.h>

const char kWindowTitle[] = "GC2A_06_ジョ_カエイ";

struct Vector3
{
	float x, y, z;
};

//加算
Vector3 Add(const Vector3& v1, const Vector3& v2) {

	return{ v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
};

//減算
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {

	return{ v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };

};

//スカラー倍
Vector3 Multiply(float& scalar, const Vector3& v) {

	return{ scalar * v.x, scalar * v.y, scalar * v.z };
};

//内積
float Dot(const Vector3& v1, const Vector3& v2) {

	return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);

};


//長さ
float Length(const Vector3& v) {

	return sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
};

//正規化
Vector3 Normalize(const Vector3& v) {
	float length = Length(v);

	if (length != 0) {
		return { v.x / length, v.y / length, v.z / length };
	}
	else {
		return { 0, 0, 0 };
	};
}


static const int kColumnWidth = 60;

void VectorScreenPrintf(int x, int y, const Vector3& v, const char* label) {
	Novice::ScreenPrintf(x, y, "%.02f", v.x);
	Novice::ScreenPrintf(x + kColumnWidth, y, "%.02f", v.y);
	Novice::ScreenPrintf(x + kColumnWidth * 2, y, "%.02f", v.z);
	Novice::ScreenPrintf(x + kColumnWidth * 3, y, "%s", label);
};

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {
	Vector3 v1{ 1.0f,3.0f,-5.0f };
	Vector3 v2{ 4.0f,-1.0f,2.0f };
	float k = { 4.0f };
	int kRowHeight = 20;


	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

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
		Vector3 resultAdd = Add(v1, v2);
		Vector3 resultSubtract = Subtract(v1, v2);
		Vector3 resultMultiply = Multiply(k, v1);
		float resultDot = Dot(v1, v2);
		float resultLength = Length(v1);
		Vector3 resultNormalize = Normalize(v2);


		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		VectorScreenPrintf(0, 0, resultAdd, " : Add");
		VectorScreenPrintf(0, kRowHeight, resultSubtract, " : Subtract");
		VectorScreenPrintf(0, kRowHeight * 2, resultMultiply, " : Multiply");
		Novice::ScreenPrintf(0, kRowHeight * 3, "%.02f : Dot", resultDot);
		Novice::ScreenPrintf(0, kRowHeight * 4, "%.02f : Length", resultLength);
		VectorScreenPrintf(0, kRowHeight * 5, resultNormalize, " : Normalize");


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
