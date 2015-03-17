vec3 ik_2_bone_centered_solver(vec3 goal, float L1, float L2)
{
#if 0 // from https://www.shadertoy.com/view/ldlGR7
	vec3 q = goal*(0.5 + 0.5*(L1*L1 - L2*L2) / dot(goal, goal));

	float s = L1*L1 - dot(q, q);
	s = max(s, 0.0);
	q += sqrt(s)*normalize(cross(goal, vec3(0, 0, 1)));

	return q;
#else // naive version with law of cosines
	float G = length(goal);

	// tetha is the angle between bone1 and goal direction
	// get it from law of cosines applied to the
	// triangle with sides: bone1, bone2, pivot_of_bone1<->goal
	float cos_theta = (L1*L1 + G*G - L2*L2) / (2.*L1*G);

	// sin^2 + cos^2 = 1 (Pythagoras in unit circle)
	float sin_theta = sqrt(1. - cos_theta * cos_theta);

	// rotation matrix by theta amount around the axis
	// perpendicular to the plane created by bone1 and bone2
	mat3 rot = mat3(
		cos_theta, -sin_theta, 0,
		sin_theta, cos_theta, 0,
		0, 0, 1.
		);

	// get the end of bone1 aka the pivot of bone2
	// by getting a vector from the goal direction
	// and rotating along with the newly found theta angle
	return rot * (normalize(goal) * L1);
#endif
}

vec3 ik_solver(vec3 start, vec3 goal, float bone_length_1, float bone_length_2)
{
	return start + ik_2_bone_centered_solver(
		goal - start, bone_length_1, bone_length_2);
}
