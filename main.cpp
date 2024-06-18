#include <Novice.h>
#include <MathFunction.h>
#include <imgui.h>



const char kWindowTitle[] = "GC2A_06_ジョ_カエイ";

bool SphereCollision(const Sphere& sphere1, const Sphere& sphere2) {
	float distance = Length(Subtract(sphere1.center, sphere2.center));
	if (distance > sphere1.radius + sphere2.radius) {
		return false;
	}
	return true;
}

static const int kRowHeight = 20;
static const int kColumnWidth = 60;


void VectorScreenPrintf(int x, int y, const Vector3& vector3, const char* name) {
	Novice::ScreenPrintf(x, y, "%6.02f  %6.02f  %6.02f  %s", vector3.x, vector3.y, vector3.z, name);
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	Vector3 rotate{ 0,0,0 };
	Vector3 translate{ 0,0,0 };
	Vector3 cameraTranslate{ 0,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };


	Sphere sphere[2];
	sphere[0].center = { 0.0f,0.0f,0.0f };
	sphere[0].radius = 0.5f;

	sphere[1].center = { 1.0f,0.0f,1.0f };
	sphere[1].radius = 0.2f;

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

		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, 1280.0f / 720.0f, 0.1f, 100.0f);
		Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
		Matrix4x4 viewPortMatrix = MakeViewportMatrix(0, 0, 1280.0f, 720.0f, 0.0f, 1.0f);

		ImGui::Begin("Window");
		ImGui::DragFloat3("sphere[0] Position", &sphere[0].center.x, 0.01f);
		ImGui::DragFloat("sphere[0] radius", &sphere[0].radius, 0.01f);
		ImGui::DragFloat3("sphere[1] Position", &sphere[1].center.x, 0.01f);
		ImGui::DragFloat("sphere[1] radius", &sphere[1].radius, 0.01f);
		ImGui::End();

	
		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///



		DrawGrid(worldViewProjectionMatrix, viewPortMatrix);
		if (SphereCollision(sphere[0], sphere[1])) {
			DrawSphere(sphere[0], worldViewProjectionMatrix, viewPortMatrix, RED);
		}
		else {
			DrawSphere(sphere[0], worldViewProjectionMatrix, viewPortMatrix, WHITE);
		}
		DrawSphere(sphere[1], worldViewProjectionMatrix, viewPortMatrix, WHITE);

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
