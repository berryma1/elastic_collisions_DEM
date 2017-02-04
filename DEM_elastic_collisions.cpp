// DEM_elastic_collisions.cpp : Defines the entry point for the console application.
// This program simulates 2D elastic collisions of spheres with one another and their container
// Created by Michael A Berry on 2/4/2017 using Visual Studio 2017

// #include "stdafx.h" // uncomment this for use in visual studio
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <time.h>
using namespace std;

// dimensions
const double step_size = 0.02;
const double PI = 3.14159;
const double sphere_diameter = 0.5;
const double sphere_mass = 1.0;

void detect_wall_collisions(vector<vector<double>> &states, int n) {
	// box height = n + 1;
	// box width = n + 1;
	for (int i = 0; i < n; i++) {
		// check left wall, then the right
		if (states[i][0] - sphere_diameter / 2 < 0) {
			// the sphere collided with the left wall
			states[i][2] *= -1;
		}
		else if (states[i][0] + sphere_diameter / 2 > n + 1) {
			// the sphere collided with the right wall
			states[i][2] *= -1;
		}

		// check bottom wall, then the top
		if (states[i][1] - sphere_diameter / 2 < 0) {
			// the sphere collided with the left wall
			states[i][3] *= -1;
		}
		else if (states[i][1] + sphere_diameter / 2 > n + 1) {
			// the sphere collided with the right wall
			states[i][3] *= -1;
		}
	}
}
void calculate_new_velocities(vector<vector<double>> &states, int i, int j) {
	cout << "Spheres " << i << " and " << j << " collided\n";
	// calculate the velocities (magnitude and direction) of spheres i and j pre-collision
	double sphere_i_vel_mag = sqrt(pow(states[i][2], 2) + pow(states[i][3], 2));
	double sphere_j_vel_mag = sqrt(pow(states[j][2], 2) + pow(states[j][3], 2));
	double sphere_i_vel_angle = atan2(states[i][3], states[i][2]); // radians
	double sphere_j_vel_angle = atan2(states[j][3], states[j][2]); // radians

	// check the kinetic energy of the spheres before and after collision to ensure conservation of energy
	double pre_collision_energy = 0.5*sphere_mass*(pow(sphere_i_vel_mag, 2) + pow(sphere_j_vel_mag, 2));

	// the two spheres that collided will exchange their velocities along the line of impact
	// calculate the collision angle: tan(theta) = delta_y / delta_x
	double collision_angle = atan2(states[i][1] - states[j][1], states[i][0] - states[j][0]); // radians

	// calculate vel_x for sphere i
	states[i][2] = states[i][2] - sphere_i_vel_mag*cos(sphere_i_vel_angle - collision_angle)*cos(collision_angle)
		+ sphere_j_vel_mag*cos(sphere_j_vel_angle - collision_angle)*cos(collision_angle);
	// calculate vel_x for sphere j
	states[j][2] = states[j][2] - sphere_j_vel_mag*cos(sphere_j_vel_angle - collision_angle)*cos(collision_angle)
		+ sphere_i_vel_mag*cos(sphere_i_vel_angle - collision_angle)*cos(collision_angle);
	// calculate vel_y for sphere i
	states[i][3] = states[i][3] - sphere_i_vel_mag*cos(sphere_i_vel_angle - collision_angle)*sin(collision_angle)
		+ sphere_j_vel_mag*cos(sphere_j_vel_angle - collision_angle)*sin(collision_angle);
	// calculate vel_y for sphere j
	states[j][3] = states[j][3] - sphere_j_vel_mag*cos(sphere_j_vel_angle - collision_angle)*sin(collision_angle)
		+ sphere_i_vel_mag*cos(sphere_i_vel_angle - collision_angle)*sin(collision_angle);

	double post_collision_energy = 0.5*sphere_mass*(pow(states[i][2], 2) + pow(states[i][3], 2) + pow(states[j][2], 2) + pow(states[j][3], 2));
}
void detect_sphere_collisions(vector<vector<double>> &states, int n, int &count) {
	double distance;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			// calculate the distance of each sphere from one another
			distance = sqrt(pow((states[i][0] - states[j][0]), 2) + pow((states[i][1] - states[j][1]), 2));
			if (distance < sphere_diameter) {
				count++;
				calculate_new_velocities(states, i, j);
			}
		}
	}
}
void write_positions(vector<vector<double>> states, int n, ofstream& file_name) {
	for (int i = 0; i < n; i++) {
		file_name << states[i][0] << ',' << states[i][1] << ',';
	}
	file_name << '\n';
}

int main()
{
	// specify the number of spheres
	int n;
	cout << "This program simulates 2D elastic collisions of spheres with one another and their container\n";
	cout << "How many spheres would you like (2-10)? ";
	cin >> n;
	if (n > 10) {
		cout << "Setting the number of spheres to the maximum: 10\n";
		n = 10;
	}
	else if (n < 2) {
		cout << "Setting the number of spheres to the minimum: 2\n";
		n = 2;
	}
	else {
		cout << "Setting the number of spheres to your selection: " << n << endl;
	}

	// represent the position and velocity of n spheres
	// sphere_states = {px, py, vx, vy}
	vector<vector<double>> sphere_states;
	sphere_states.resize(n, vector<double>(4, 0));
	
	for (int i = 0; i < n; i++) {
		// arrange the spheres diagonally in the container
		sphere_states[i][0] = double(i + 1);
		sphere_states[i][1] = double(i + 1);
		// randomly generate the direction of initial velocity vector for each sphere
		srand(time(NULL));
		double temp = double(rand() % 360); // angle in degrees
		sphere_states[i][2] = cos(temp * PI / 180.0);
		sphere_states[i][3] = sin(temp * PI / 180.0);
	}
	// create a file to write to
	ofstream positions;
	positions.open("elastic_collisions_res.csv");

	// setup complete
	int number_of_collisions = 0; // sphere to sphere collisions
	int number_of_steps = 0;
	while (number_of_collisions < 10) {
		// calculate the new position of each particle at the end of each time step
		for (int i = 0; i < n; i++) {
			sphere_states[i][0] += step_size*sphere_states[i][2];
			sphere_states[i][1] += step_size*sphere_states[i][3];
		}
		detect_wall_collisions(sphere_states, n);
		detect_sphere_collisions(sphere_states, n, number_of_collisions);
		write_positions(sphere_states, n, positions);
		number_of_steps++;
	}
	return 0;
}
