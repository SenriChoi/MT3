#define NOMINMAX
#include <Novice.h>
#include <MathFunction.h>
#include <imgui.h>
#include <algorithm>

const char kWindowTitle[] = "GC2A_06_ジョ_カエイ";



struct AABB {
	Vector3 min;
	Vector3 max;
};

bool IsCollision(const AABB& aabb, const Sphere& sphere) {
	Vector3 clossestPoint
	{
		std::clamp(sphere.center.x,aabb.min.x,aabb.max.x),
		std::clamp(sphere.center.y,aabb.min.y,aabb.max.y),
		std::clamp(sphere.center.z,aabb.min.z,aabb.max.z)
	};
	
	float distance = Length(Subtract(clossestPoint, sphere.center));

	if (distance <= sphere.radius) {
		return true;
	}

	return false;
}

void CorrectAABB(AABB& aabb) {
	aabb.max.x = std::max(aabb.max.x, aabb.min.x);
	aabb.min.x = std::min(aabb.min.x, aabb.max.x);
	aabb.max.y = std::max(aabb.max.y, aabb.min.y);
	aabb.min.y = std::min(aabb.min.y, aabb.max.y);
	aabb.max.z = std::max(aabb.max.z, aabb.min.z);
	aabb.min.z = std::min(aabb.min.z, aabb.max.z);
}


void DrawAABB(const AABB& aabb, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewPortMatrix, uint32_t color) {

	Vector3 AABBpoint[8];
	AABBpoint[0] = Vector3(aabb.min.x, aabb.min.y, aabb.min.z);
	AABBpoint[1] = Vector3(aabb.max.x, aabb.min.y, aabb.min.z);
	AABBpoint[2] = Vector3(aabb.max.x, aabb.min.y, aabb.max.z);
	AABBpoint[3] = Vector3(aabb.min.x, aabb.min.y, aabb.max.z);
	AABBpoint[4] = Vector3(aabb.min.x, aabb.max.y, aabb.min.z);
	AABBpoint[5] = Vector3(aabb.max.x, aabb.max.y, aabb.min.z);
	AABBpoint[6] = Vector3(aabb.max.x, aabb.max.y, aabb.max.z);
	AABBpoint[7] = Vector3(aabb.min.x, aabb.max.y, aabb.max.z);

	Vector3 screenPoints[8];
	for (uint32_t index = 0; index < 8; ++index) {
		Vector3 ndcVertex = Transform(AABBpoint[index], viewProjectionMatrix);
		screenPoints[index] = Transform(ndcVertex, viewPortMatrix);
	}

	Novice::DrawLine((int)screenPoints[0].x, (int)screenPoints[0].y, (int)screenPoints[1].x, (int)screenPoints[1].y, color);
	Novice::DrawLine((int)screenPoints[1].x, (int)screenPoints[1].y, (int)screenPoints[2].x, (int)screenPoints[2].y, color);
	Novice::DrawLine((int)screenPoints[2].x, (int)screenPoints[2].y, (int)screenPoints[3].x, (int)screenPoints[3].y, color);
	Novice::DrawLine((int)screenPoints[3].x, (int)screenPoints[3].y, (int)screenPoints[0].x, (int)screenPoints[0].y, color);
	Novice::DrawLine((int)screenPoints[4].x, (int)screenPoints[4].y, (int)screenPoints[5].x, (int)screenPoints[5].y, color);
	Novice::DrawLine((int)screenPoints[5].x, (int)screenPoints[5].y, (int)screenPoints[6].x, (int)screenPoints[6].y, color);
	Novice::DrawLine((int)screenPoints[6].x, (int)screenPoints[6].y, (int)screenPoints[7].x, (int)screenPoints[7].y, color);
	Novice::DrawLine((int)screenPoints[7].x, (int)screenPoints[7].y, (int)screenPoints[4].x, (int)screenPoints[4].y, color);
	Novice::DrawLine((int)screenPoints[0].x, (int)screenPoints[0].y, (int)screenPoints[4].x, (int)screenPoints[4].y, color);
	Novice::DrawLine((int)screenPoints[1].x, (int)screenPoints[1].y, (int)screenPoints[5].x, (int)screenPoints[5].y, color);
	Novice::DrawLine((int)screenPoints[2].x, (int)screenPoints[2].y, (int)screenPoints[6].x, (int)screenPoints[6].y, color);
	Novice::DrawLine((int)screenPoints[3].x, (int)screenPoints[3].y, (int)screenPoints[7].x, (int)screenPoints[7].y, color);

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

	AABB aabb{ {-0.5f,-0.5f,-0.5f},{0.0f,0.0f,0.0f} };
	Sphere sphere{ {0.0f,0.0f,0.0f},0.5f };


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

		CorrectAABB(aabb);


		ImGui::Begin("Window");
		ImGui::DragFloat3("Sphere center", &sphere.center.x, 0.01f);
		ImGui::DragFloat("Sphere radius", &sphere.radius, 0.01f);
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::End();

	
		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		/// 
		DrawGrid(worldViewProjectionMatrix, viewPortMatrix);


		DrawSphere(sphere, worldViewProjectionMatrix, viewPortMatrix, WHITE);

		if (IsCollision(aabb, sphere)) {
			DrawAABB(aabb, worldViewProjectionMatrix, viewPortMatrix, RED);
		}
		else {
			DrawAABB(aabb, worldViewProjectionMatrix, viewPortMatrix, WHITE);
		}
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
