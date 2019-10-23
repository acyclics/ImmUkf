#include "mujoco.h"
#include "glfw3.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <chrono>
#include <iostream>

#include "../armorfilter.h"
#include "../utils.h"
#include "../ukf_plot.h"


// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
mjvCamera cam;                      // abstract camera
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     // custom GPU context

// mouse interaction
bool button_left = false;
bool button_middle = false;
bool button_right = false;
double lastx = 0;
double lasty = 0;

// keyboard callback
void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods)
{
	// backspace: reset simulation
	if (act == GLFW_PRESS && key == GLFW_KEY_BACKSPACE)
	{
		mj_resetData(m, d);
		mj_forward(m, d);
	}
}

// mouse button callback
void mouse_button(GLFWwindow* window, int button, int act, int mods)
{
	// update button state
	button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS);
	button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS);
	button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS);

	// update mouse position
	glfwGetCursorPos(window, &lastx, &lasty);
}


// mouse move callback
void mouse_move(GLFWwindow * window, double xpos, double ypos)
{
	// no buttons down: nothing to do
	if (!button_left && !button_middle && !button_right)
		return;

	// compute mouse displacement, save
	double dx = xpos - lastx;
	double dy = ypos - lasty;
	lastx = xpos;
	lasty = ypos;

	// get current window size
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// get shift key state
	bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
		glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

	// determine action based on mouse button
	mjtMouse action;
	if (button_right)
		action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
	else if (button_left)
		action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
	else
		action = mjMOUSE_ZOOM;

	// move camera
	mjv_moveCamera(m, action, dx / height, dy / height, &scn, &cam);
}


// scroll callback
void scroll(GLFWwindow * window, double xoffset, double yoffset)
{
	// emulate vertical mouse motion = 5% of window height
	mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05 * yoffset, &scn, &cam);
}


// main function
int main()
{
	// activate software
	mj_activate("C:/Users/impec/Desktop/Mods/ai/mujoco/mujoco200_win64/bin/mjkey.txt");

	// load and compile model
	char error[1000] = "Could not load binary model";
	m = mj_loadXML("./mujoco/armor.xml", 0, error, 1000);

	// make data
	d = mj_makeData(m);

	// init GLFW
	if (!glfwInit())
		mju_error("Could not initialize GLFW");

	// create window, make OpenGL context current, request v-sync
	GLFWwindow* window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	// initialize visualization data structures
	mjv_defaultCamera(&cam);
	mjv_defaultOption(&opt);
	mjv_defaultScene(&scn);
	mjr_defaultContext(&con);

	// create scene and context
	mjv_makeScene(m, &scn, 2000);
	mjr_makeContext(m, &con, mjFONTSCALE_150);

	// install GLFW mouse and keyboard callbacks
	glfwSetKeyCallback(window, keyboard);
	glfwSetCursorPosCallback(window, mouse_move);
	glfwSetMouseButtonCallback(window, mouse_button);
	glfwSetScrollCallback(window, scroll);
	
	// plot
	const int plot_size = 100;
	vector<double> residue(plot_size), time(plot_size), var(plot_size);
	int int_t = 0;

	armorfilter af;
	double spd = 3;
	double predicted_dist2center = 0.5;
	Vector3d new_pos;

	cout << "Begin\n";

	// run main loop, target real-time simulation and 60 fps rendering
	while (!glfwWindowShouldClose(window))
	{
		// advance interactive simulation for 1/60 sec
		//  Assuming MuJoCo can simulate faster than real-time, which it usually can,
		//  this loop will finish on time for the next frame to be rendered at 60 fps.
		//  Otherwise add a cpu timer and exit this loop when it is time to render.
		mjtNum simstart = d->time;
		while (d->time - simstart < 1.0 / 60.0) {
			mj_step(m, d);
		}

		// get framebuffer viewport
		mjrRect viewport = { 0, 0, 0, 0 };
		glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

		// update scene and render
		mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
		mjr_render(viewport, &scn, &con);
		
		// Control armor
		d->ctrl[0] = 7;
		//d->ctrl[1] = spd;

		if (int_t % 100 == 0) spd = -spd;

		// Hold indicator at height
		d->qfrc_applied[6] = d->qfrc_bias[6];

		Eigen::Quaterniond q(d->xquat[16], d->xquat[17], d->xquat[18], d->xquat[19]);
		q = q.normalized();
		auto rollpitchyaw = q.toRotationMatrix().eulerAngles(0, 1, 2);

		if (int_t % 10 == 0) {
			sensordata sd1 = { d->xpos[15] , d->xpos[16], d->xpos[17], d->xquat[16], d->xquat[17], d->xquat[18], d->xquat[19] };
			sensordata sd2 = { d->xpos[18] , d->xpos[19], d->xpos[20], d->xquat[20], d->xquat[21], d->xquat[22], d->xquat[23] };

			vector<sensordata> vsd = { sd1, sd2 };
			double time = af.get_tda().current_time();
			//vector<double> vtime = { time };
			vector<double> vtime = { (double)int_t, (double)int_t };
			af.update(vsd, vtime);
		}
		else {
			//af.predict();
			af.predict(int_t);
		}

		VectorXd X = af.get_tda().get_trackers()[0].get_imm().get();
		//VectorXd X = af.get_tda().get_trackers()[0].get_imm().peek(15);
		//VectorXd X = af.get_tda().get_trackers()[0].get_imm().get_xi(0);

		VectorXd mu = af.get_tda().get_trackers()[0].get_imm().get_mu();
		for (int i = 0; i < NM; ++i) {
			cout << "Model " << i << ": " << mu(i) << "    |    ";
		}
		cout << "\n";
		
		cout << "x: " << X[0] << "   Actual: " << d->xpos[15] << "\n"
			<< "y: " << X[1] << "   Actual: " << d->xpos[16] << "\n"
			<< "z: " << X[2] << "   Actual: " << d->xpos[17] << "\n"
			<< "roll: " << X[3] << "   Actual: " << rollpitchyaw[0] << "\n"
			<< "pitch: " << X[4] << "   Actual: " << rollpitchyaw[1] << "\n"
			<< "yaw: " << X[5] << "   Actual: " << rollpitchyaw[2] << "\n";

		cout << "\n\n\n";

		d->qpos[4] = X[0];
		d->qpos[5] = X[1];
		d->qpos[6] = X[2];

		/*
		residue[int_t] = d->xpos[15] - imm.get_ukf(1).get_measurement_pred(0);
		time[int_t] = int_t;
		var[int_t] = imm.get_ukf(1).get_state_var(0, 0);// +imm.get_ukf(1).get_noise(0);

		if (plot_size - 1 == int_t) {
			write_residueVsTime(time, residue, var, plot_size);
			while (1);
		}
		*/

		// swap OpenGL buffers (blocking call due to v-sync)
		glfwSwapBuffers(window);

		// process pending GUI events, call GLFW callbacks
		glfwPollEvents();
		
		++int_t;
	}

	//free visualization storage
	mjv_freeScene(&scn);
	mjr_freeContext(&con);

	// free MuJoCo model and data, deactivate
	mj_deleteData(d);
	mj_deleteModel(m);
	mj_deactivate();

	// terminate GLFW (crashes with Linux NVidia drivers)
#if defined(__APPLE__) || defined(_WIN32)
	glfwTerminate();
#endif

	return 1;
}
