cbuffer uniforms : register(b0)
{
	float2 u_res;
	float u_time;
	float2 u_mouse;
};

float4 PS_main(float4 uv : SV_Position) : SV_Target
{
	return float4(1, 0, 1, 1);
}