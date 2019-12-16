#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

/*
 * File:          nav0.c
 * Date:   1 Tir 1398
 * Author: Safa Mohammadi, parisa ghasemi, mohaddese Fahiminia
 */


#include <webots/differential_wheels.h>
#include <webots/distance_sensor.h>
#include <webots/robot.h>
#include <webots/motor.h>
#include <webots/led.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <ctype.h>
#include <assert.h>

typedef char** Map;

#define PI 3.14159
Map primaryMap;
double mapXScale = 0.5; 
double mapYScale = 0.5; 
const int mapHeight = 192; 
const int mapWidth = 174; 
const int robot_angle =3* PI / 2.0;
const int numOfParticles = 3000;
double Xmax = mapWidth * mapXScale;
double Ymax = mapHeight * mapYScale;


#define WHEEL_RADIUS 0.02
#define AXLE_LENGTH 0.052
#define TIME_STEP 128 
#define SPEED 1.0 
#define ROTATION_SPEED 0.5


double sensor_values[8];
typedef struct { double lms; double rms; int tss; double mean; double variance; } Action; //left motor speed,right motor speed, time steps
Action* Actions;
/* motion model
distances: 2		5		10
means:	   2.27		5.1     9.9
variances: 0.0161   0.11    0.11
*/

#define syncadd "sync.txt"
#define datafileadd "datafile.txt"

typedef struct { float x; float y; float theta; } state;

double distance_means[12] = { 0, 1, 2, 3, 4, 5.1, 6.2, 7.2, 8.2, 9.3, 10.30, 20.60 };
double distance_variances[12] = { .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5 };
double orientation_mean = 0;
double orientation_var = 0.017;

double sensor_means[8][3] = { 
{ 1568.9, 581.1, 510.8},
{ 2300.7, 333.7, 244.2},
{ 3325.6, 497.7, 497.7},
{ 1621.06, 941.5, 878.8},
{ 1629.6, 781.7, 721.07},
{ 3285.3, 490.6, 397.6},
{ 2261.6, 359.7, 267.2},
{ 1750.09, 711.6, 634.5} };

double sensor_variances[8][3] = { 
{ 17.4, 10.4, 11.08},
{ 8.2, 3.5, 7.5},
{ 2.3, 0.8, 0.8},
{ 12.9, 5.6, 5.4},
{ 1.5, 2.3, 5.3},
{ 122.9, 1.8, 3.3},
{ 1.6, 3.2, 5.6},
{ 23.5, 6.6, 7.8} };

const double sigmaPosision[] = { 0.3, 0.3, 0.01 };
double wSlow = 0.002;
double wFast = 0.006;
double alphaSlow = 0.1;
double alphaFast = 0.3;

typedef struct
{
	double x;
	double y;
	double theta;
	double weight;
}Particle;

Particle** allParticles;

double sensor_nearest_obstacle(Particle* particle, int index);
double estimate_sensor_mean(double dist, int sindex);
double estimate_sensor_var(double dist, int sindex);
float pdf_gaussian(double x, double mu, double var);
bool check_for_convergence();
double sample_gaussian(double mean, double variance);
bool does_intersect(double x, double y);

void show(int rounds, FILE* sync, FILE* datafile) {
	if (datafile == NULL) {
		printf("failure open datafile in show , %d\n", rounds);
	}
	else {
		int i;
		for (i = 0; i < numOfParticles; i++) {
			fprintf(datafile, "%f %f %f\n", floor(allParticles[i]->x / mapXScale),
				floor(allParticles[i]->y / mapYScale), allParticles[i]->theta);
		}
		fprintf(datafile, "END\n");
		fflush(datafile);
	}
	if (sync == NULL) {
		printf("failure open sync in show , %d\n", rounds);
	}
	else {
		printf("rounds = %d\n", rounds);
		fprintf(sync, "%d\n", rounds);
		fflush(sync);
	}
}

int isThereColision(double x, double y) {
	int xMap = (int)(y / mapXScale);
	int ymap = (int)(x / mapYScale);
	if (xMap < 3 || xMap >(numOfParticles - 3) || ymap < 3 || ymap >(numOfParticles - 3))
		return 1;
	if (
		primaryMap[xMap - 1][ymap - 1] == '0' &&
		primaryMap[xMap][ymap - 1] == '0' &&
		primaryMap[xMap + 1][ymap - 1] == '0' &&
		primaryMap[xMap - 1][ymap] == '0' &&
		primaryMap[xMap][ymap] == '0' &&
		primaryMap[xMap + 1][ymap] == '0' &&
		primaryMap[xMap - 1][ymap + 1] == '0' &&
		primaryMap[xMap][ymap + 1] == '0' &&
		primaryMap[xMap + 1][ymap + 1] == '0')
		return 0;
	else
		return 1;
}

Particle * generateParticle(int n, int check) {
	Particle* newParticle = (Particle*)malloc(sizeof(Particle));

	while (1){
		newParticle->x = (rand() % ((int)Xmax)) + sample_gaussian(0.0, 1.0);
		newParticle->y = (rand() % ((int)Ymax)) + sample_gaussian(0.0, 1.0);
		newParticle->theta = (rand() % 360) * PI / 180;
		if (newParticle->x <= Xmax || newParticle->x > 0 || newParticle->y > 0 || newParticle->y <= Ymax) {
			if (!(does_intersect(newParticle->x, newParticle->y))) {
				break;
			}
		}
	}
	if (n == 1)
		newParticle->weight = (1000.0 / (rand() % (numOfParticles)+1));
	else
		newParticle->weight = (1000.0 / (n));
	newParticle->weight *= 0.001;

	return newParticle;
}

void uniformSampling(int n) {
	Particle* p;
	int check = 0;
	for (int i = 0; i < n; i++) {
		while (1) {
			p = generateParticle(n, check);
			check = 1 - check;
			free(p);
		}
		allParticles[i] = p;
	}
}

Particle** newParticles;
FILE* outwe;
bool wfset = false;

void resampling() {
	double avgWeight = 0;
	double maxWeight = 0;
	for (int i = 0; i < numOfParticles; i++) {
		avgWeight += allParticles[i]->weight / numOfParticles;
		if (allParticles[i]->weight > maxWeight) {
			maxWeight = allParticles[i]->weight;
		}
	}
	if (!wfset) {
		outwe = fopen("ooo.txt", "w+");
	}
	double sum = 0;
	for (int i = 0; i < numOfParticles; i++)
		sum += allParticles[i]->weight;
	for (int i = 0; i < numOfParticles; i++)
	{
		allParticles[i]->weight = allParticles[i]->weight / sum;
		fprintf(outwe, "%f", (float)allParticles[i]->weight);
	}
	fprintf(outwe, "\n");

	wSlow += alphaSlow * (avgWeight - wSlow);
	wFast += alphaFast * (avgWeight - wFast);

	double beta = 0.0;
	int particlesSize = sizeof(allParticles) / sizeof(Particle*);

	double random_threshold = RAND_MAX * (1 - wFast / wSlow);
	int index = 0;
	double eps = ((double)rand()) / ((double)RAND_MAX);
	double unit = 1.0 / ((double)numOfParticles);
	double sw = 0.0;
	Particle * p;
    for (int i = 0; i < numOfParticles - 20; i++)
	{
		int z = rand();
        while (sw < eps) {
            index = (index + 1) % numOfParticles;
            sw += allParticles[index]->weight;
        }

        newParticles[i]->x = allParticles[index]->x;
        newParticles[i]->y = allParticles[index]->y;
        newParticles[i]->theta = allParticles[index]->theta;
        newParticles[i]->weight = allParticles[index]->weight; //!!!!!!
		eps += unit;
	}

	for (int i = numOfParticles - 20; i < numOfParticles; i++)
	{
		newParticles[i] = generateParticle(2 * numOfParticles, 0);
	}
	double sumOfWeights = 0.0;
	for (int i = 0; i < numOfParticles; i++) {
		sumOfWeights += newParticles[i]->weight;
	}
	for (int i = 0; i < numOfParticles; i++) {
		newParticles[i]->weight /= sumOfWeights;
	}
	Particle** temp = allParticles;
	allParticles = newParticles;
	newParticles = allParticles;
}

double sample_gaussian(double mean, double variance) {
	double r = 0;
	int v = (int)(sqrtl(variance) * 1000);
	for (int i = 0; i < 12; i++) {
		r += ((rand() % (2 * v)) - v) / 1000.0;
	}
	return mean + (r / 2.0);
}

bool does_intersect(double x, double y) {
	double r = 2;
	double bottomleftx = (int)((x - r) / mapXScale);
	double bottomlefty = (int)((y - r) / mapYScale);
	double toprightx = (int)((x + r) / mapXScale);
	double toprighty = (int)((y + r) / mapYScale);

	for (int i = bottomleftx; i <= toprightx; i++) {
		for (int j = bottomlefty; j <= toprighty; j++) {
			if (i >= 0 && j >= 0 && i < mapWidth && j < mapHeight && primaryMap[j][i] == '1') {
				return true;
			}
		}
	}
	return false;
}

bool path_intersects(double sx, double sy, double dx, double dy, Map primaryMap) {
	if (does_intersect(sx, sy) || does_intersect(sx + dx, sy + dy)) {
		return true;
	}
	int numpoints = 1;
	double x = sx, y = sy;
	for (int i = 0; i < numpoints; i++) {
		x += (i + 1) * dx / (numpoints + 1.0);
		y += (i + 1) * dy / (numpoints + 1.0);
		if (does_intersect(x, y)) {
			return true;
		}
	}
	return false;
}

void motion_model(Particle * p, Action a) {
	if (a.lms == a.rms) { 
		double d = sample_gaussian(a.mean, a.variance);
		double angle = p->theta + sample_gaussian(0, 0.0001);
		double dx = (-1 * d * sin(angle));
		double dy = (d * cos(angle));
		if (path_intersects(p->x, p->y, dx, dy, primaryMap)) {
			p->weight = 0; 
		}
		p->x += dx;
		p->y += dy;
		p->theta = angle;
	}
	else {
		double angle = 0;
		if (a.lms > a.rms) { 
			angle = (-1 * (PI / 2)) + sample_gaussian(0, 0.0001);
			p->theta += angle;
		}
		else if (a.lms < a.rms) {
			angle = (PI / 2) + sample_gaussian(0, 0.0001);
			p->theta += angle;

		}
	}
	if (p->x >= Xmax || p->x < 0 || p->y < 0 || p->y >= Ymax) {
		p->weight = 0;
	}
	if (does_intersect(p->x, p->y)) {
		p->weight = 0.0;\
	}
}

void apply_motion_model(Action a) {
	for (int i = 0; i < numOfParticles; i++) {
		motion_model(allParticles[i], a);
	}
}

double estimate_sensor_mean(double dist, int sindex) { 
	double mean, var;
	double a, b;
	if (dist < 1) {
		a = ((sensor_means[sindex][4] - sensor_means[sindex][1]) * 2 / 3.0) - ((sensor_means[sindex][7] - sensor_means[sindex][4]) / 3.0);
		b = sensor_means[sindex][1] - (a * 1.0);
	}
	else if (dist>=1 && dist < 4) {
		a = ((sensor_means[sindex][4] - sensor_means[sindex][1]) / 3.0);
		b = sensor_means[sindex][4] - (a * 4.0);
	}
	else if (dist < 7) {
		a = ((sensor_means[sindex][7] - sensor_means[sindex][4]) / 3.0);
		b = sensor_means[sindex][4] - (a * 4.0);
	}
	else {
		a = ((sensor_means[sindex][7] - sensor_means[sindex][4]) * 2 / 3.0) - ((sensor_means[sindex][4] - sensor_means[sindex][1]) / 3.0);;
		b = sensor_means[sindex][7] - (a * 4.0);
	}
	return (a * dist) - b;
}

double estimate_sensor_var(double dist, int sindex) { 
	double mean, var;
	double a, b;
	if (dist < 1) {
		a = ((sensor_variances[sindex][4] - sensor_variances[sindex][1]) * 2 / 3.0) - ((sensor_variances[sindex][7] - sensor_variances[sindex][4]) / 3.0);
		b = sensor_variances[sindex][1] - (a * 1.0);
	}
	else if (dist < 4) {
		a = ((sensor_variances[sindex][4] - sensor_variances[sindex][1]) / 3.0);
		b = sensor_variances[sindex][4] - (a * 4.0);
	}
	else if (dist < 7) {
		a = ((sensor_variances[sindex][7] - sensor_variances[sindex][4]) / 3.0);
		b = sensor_variances[sindex][4] - (a * 4.0);
	}
	else {
		a = ((sensor_variances[sindex][7] - sensor_variances[sindex][4]) * 2 / 3.0) - ((sensor_variances[sindex][4] - sensor_variances[sindex][1]) / 3.0);;
		b = sensor_variances[sindex][7] - (a * 4.0);
	}
	return (a * dist) - b;
}

double sensorModel(Particle* p) { 
	double prob = 100000000;
	for (int i = 0; i < 8; i++) {
		double d = sensor_nearest_obstacle(p, i);
		if (d < 1 || d>7)
			continue;
		double mu = estimate_sensor_mean(d, i);
		double v = estimate_sensor_var(d, i);
		prob *= pdf_gaussian(sensor_values[i],mu, v);
	}
	return prob;
}

void apply_sensor_model() {
	double weight_sum = 0.0;
    for (int i = 0; i < numOfParticles; i++) {
		double p = sensorModel(allParticles[i]);
		allParticles[i]->weight = p * (allParticles[i]->weight);
		weight_sum += allParticles[i]->weight;
	}
	if (weight_sum != 0) {
		for (int i = 0; i < numOfParticles; i++) {
			allParticles[i]->weight /= weight_sum;
		}
	}
	else {
		printf("All weights equal zero!!!!!!!!!!!!!!!!!\n\n\n\n");
	}
}

void readmap() {

	primaryMap = (char**)malloc(mapHeight * sizeof(char*));
	for (int i = 0; i < mapHeight; i++) {
		primaryMap[i] = (char*)malloc(mapWidth * sizeof(char));
		for (int j = 0; j < mapWidth; j++) {
			primaryMap[i][j] = '0';
		}
	}

	FILE* f = fopen("map.txt", "r");
	FILE* fo = fopen("out.txt", "w");
	char c;
	c = (char)fgetc(f);
	int line = mapHeight - 1;
	int col = 0;
	while (c != EOF) {
		if (c == '\n' || c == '\r') {
			fprintf(fo, "\n");
			line--;
			col = 0;
		}
		if (c == '0' || c == '1')
		{
			primaryMap[line][col] = c;
			fprintf(fo, "%c ", primaryMap[line][col]);
			col++;
		}
		c = (char)fgetc(f);
	}
	fclose(fo);
	fclose(f);
}





double sensor_nearest_obstacle(Particle* particle, int index) {
	double sensor_angles[8] = { 1.27 - PI / 2, 0.77 - PI / 2, 0 - PI / 2, 5.21 - PI / 2,
		4.21 - PI / 2, 3.1415 - PI / 2, 2.37 - PI / 2, 1.87 - PI / 2 };
	double theta = particle ->theta + sensor_angles[index];
	double cot1 = cos(theta) / sin(theta);
	double sin1 = sin(theta);
	int x = (int)(particle ->x / mapXScale);
	int y = (int)(particle ->y / mapYScale);
	int i = 0, j = 0;
	while ((y + j) < mapHeight && (x + i) < mapWidth && (x + i) >= 0 && (y + j) >= 0 && primaryMap[y + j][x + i] == '0') {
		if (cot1 > 0) {
			if ((double)i / j > cot1) {
				if (sin1 > 0) { j--; }
				else { j++; }
			}
			else {
				if (sin1 > 0) { i++; }
				else { i--; }
			}
		}
		else {
			if ((float)i / j > abs(cot1)) {
				if (sin1 > 0) { j--; }
				else { j++; }
			}
			else {
				if (sin1 > 0) { i--; }
				else { i++; }
			}
		}
	}
	i *= mapXScale;
	j *= mapYScale;
	double distance = sqrt(i * i + j * j);
	return distance;
}



Action * init_action(Action * act, double lms, double rms, int tss, double mean, double variance) {
	Action* a = act;
	if (a == NULL) {
		a = (Action*)malloc(sizeof(Action));
	}
	a->lms = lms;
	a->rms = rms;
	a->tss = tss;
	a->mean = mean;
	a->variance = variance;
	return a;
}

void init_motions() { 
	Actions = (Action*)malloc(sizeof(Action) * 5);
	init_action(&Actions[0], -1.0 * SPEED, SPEED, (int)(PI * 1000 / (2.0 * SPEED * TIME_STEP)), PI / 2.0, 0.1); //LEFTTURN
	init_action(&Actions[1], SPEED, -1 * SPEED, (int)(PI * 1000 / (2.0 * SPEED * TIME_STEP)), PI / 2.0, 0.1); //RIGHTTURN
	init_action(&Actions[2], SPEED, SPEED, (int)(0.62 * 2.0 * 1000 / (SPEED * TIME_STEP)), 2.27, sqrt(0.0161)); //FORWARD2  
	init_action(&Actions[3], -1 * SPEED, -1 * SPEED, (int)(0.62 * 5.0 * 1000 / (SPEED * TIME_STEP)), 2.27, sqrt(0.0161)); //BACKWARD2  
	init_action(&Actions[3], SPEED, SPEED, (int)(0.6 * 2.0 * 1000 / (SPEED * TIME_STEP)), 5.1, sqrt(0.11)); //FORWARD5
	init_action(&Actions[4], SPEED, SPEED, (int)(0.63 * 2.0 * 1000 / (SPEED * TIME_STEP)), 9.9, sqrt(0.11)); //FORWARD10
}


state * new_state(float x, float y, float theta) {
	state* s = (state*)malloc(sizeof(state));
	s->theta = theta;
	s->x = x;
	s->y = y;
}

void printAction(Action * a) {
	printf("leftmotor speed = %f\nrightmotorspeed = %f\ntimesteps = %d\n,meanDistance = %f\ndistance variance = %f\n***",
		a->lms, a->rms, a->tss, a->mean, a->variance);
}

void copy_act(Action * a1, Action * a2) {
	init_action(a2, a1->lms, a1->rms, a1->tss, a1->mean, a1->variance);
}

void random_action(Action * act) { 
	int x = (int)(rand() % 3); 
	copy_act(&Actions[x], act);
}

void set_action(Action * act, double sensor_values[]) {
	static double max_readings[8] = { 400,500,600,532,625,600,300,480 };

	bool openWays[8];
	for (int i = 0; i < 8; i++) {
		openWays[i] = (sensor_values[i] < max_readings[i]);
	}
	if (openWays[0] && openWays[7]) { 
		copy_act(&Actions[2], act); 
	}
	else if (openWays[1]) {
		copy_act(&Actions[1], act);
	}
	else if (openWays[6]) {
		copy_act(&Actions[0], act); 
	}
	else { 
		copy_act(&Actions[3], act); 
	}
	int x = (int)(rand());
	if (x % 10 == 0) {
		random_action(act);
	}
}

float pdf_gaussian(double x, double mu, double var) {
	return 0.01* (100 / sqrt(2 * PI * var * var)) * (1 + erf(-1 * (x - mu) * (x - mu) / (2 * (var * var))));
}

void read_distance_sensors(WbDeviceTag tags[]) {
	for (int i = 0; i < 8; i++) {
		sensor_values[i] = wb_distance_sensor_get_value(tags[i]);
		printf("s%d : %f //", i, (float)sensor_values[i]);
	}
	printf("\n");
}


int main(int argc, char** argv)
{
	newParticles = (Particle * *)malloc(sizeof(Particle*) * numOfParticles);
	for (int i = 0; i < numOfParticles; i++)
	{
		Particle* p = (Particle*)malloc(sizeof(Particle));;
		newParticles[i] = p;
	}
	allParticles = (Particle * *)malloc(sizeof(Particle*) * numOfParticles);

	WbDeviceTag left_motor, right_motor;
	Action* act = (Action*)malloc(sizeof(Action));
	WbDeviceTag distance_sensor[8];
	int i;
	bool firstact = true;

	//initialization
	srand(time(0));
	init_motions();

	FILE* sync = fopen(syncadd, "w");
	FILE* datafile = fopen(datafileadd, "w");
	readmap();

	wb_robot_init();
	left_motor = wb_robot_get_device("left wheel motor");
	right_motor = wb_robot_get_device("right wheel motor");
	wb_motor_set_position(left_motor, INFINITY);
	wb_motor_set_position(right_motor, INFINITY);
	wb_motor_set_velocity(left_motor, 0.0);
	wb_motor_set_velocity(right_motor, 0.0);


	for (i = 0; i < 8; i++) {
		char device_name[4];

		/* get distance sensors */
		sprintf(device_name, "ps%d", i);
		distance_sensor[i] = wb_robot_get_device(device_name);
		wb_distance_sensor_enable(distance_sensor[i], TIME_STEP);
	}

	//get leds
	WbDeviceTag leds[8];
	leds[0] = wb_robot_get_device("led0");
	leds[1] = wb_robot_get_device("led1");
	leds[2] = wb_robot_get_device("led2");
	leds[3] = wb_robot_get_device("led3");
	leds[4] = wb_robot_get_device("led4");
	leds[5] = wb_robot_get_device("led5");
	leds[6] = wb_robot_get_device("led6");
	leds[7] = wb_robot_get_device("led7");


	int waittime = 0;
	uniformSampling(numOfParticles);
	int rounds = 1;
	show(rounds, sync, datafile);
	while (wb_robot_step(TIME_STEP) != -1) {
		if (act->tss == 0 || firstact) {
			wb_motor_set_velocity(left_motor, 0);
			wb_motor_set_velocity(right_motor, 0);
			if (waittime >0) {
				waittime--;
				bool waiting = true;
				continue;
			}
			firstact = false;
			read_distance_sensors(distance_sensor);
			if(rounds > 3) {
				apply_sensor_model();
				resampling();
			}
			if (check_for_convergence()) {
				wb_motor_set_velocity(left_motor, 0);
				wb_motor_set_velocity(right_motor, 0);

				wb_led_set(leds[0], 2);
				wb_led_set(leds[1], 2);
				wb_led_set(leds[2], 2);
				wb_led_set(leds[3], 2);
				wb_led_set(leds[4], 2);
				wb_led_set(leds[5], 2);
				wb_led_set(leds[6], 2);
				wb_led_set(leds[7], 2);
				continue;
			}
			rounds++;
			show(rounds, sync, datafile);
			if (act->lms < 0 && act->rms < 0) {  
				copy_act(&Actions[1], act); 
			}
			else {
				set_action(act, sensor_values);
			}
			waittime = 10;
			apply_motion_model(*act);
			show(rounds, sync, datafile);
			continue;
		}


		wb_motor_set_velocity(left_motor, act->lms);
		wb_motor_set_velocity(right_motor, act->rms);
		act->tss -= 1;
	}

	free(act);
	fclose(sync);
	fclose(datafile);
	fclose(outwe);
	wb_robot_cleanup();
	return 0;
}

bool check_for_convergence()
{
	bool hasConverged = false;
	int belief = 0;
	double maxWeight = 0;
	double avgx = 0;
	double avgy = 0;
	for (int i = 0; i < numOfParticles; i++)
	{
		if (allParticles[i]->weight > maxWeight) {
			maxWeight = allParticles[i]->weight;
			belief = i;
		}
		avgx += allParticles[i]->x;
		avgy = allParticles[i]->y;
	}

	avgx /= (double)numOfParticles;
	avgy /= (double)numOfParticles;
	int convergedParticles = 0;
	double validDistance = 20.00;
	double distance;
	for (int i = 0; i < numOfParticles; i++)
	{
		distance = (allParticles[i]->x - avgx) * (allParticles[i]->x - avgx) +
			(allParticles[i]->y - avgy) * (allParticles[i]->y - avgy);
		if (distance < validDistance)
			convergedParticles++;

	}
	if (convergedParticles > 0.8 * numOfParticles)
	{
		hasConverged = true;
		printf("\nThe robot is in cordinates (%f,%f);\n", allParticles[belief]->x, allParticles[belief]->y);
	}
	return hasConverged;
}
