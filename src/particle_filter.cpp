/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	num_particles = 30;
        default_random_engine gen;	
	double sigma_x,sigma_y,sigma_theta;
	sigma_x = std[0];
	sigma_y = std[1];
	sigma_theta = std[2];

	normal_distribution<double> dist_x(x, sigma_x);
	normal_distribution<double> dist_y(y, sigma_y);
	normal_distribution<double> dist_theta(theta, sigma_theta);

	for (int i=0;i<num_particles;i++) {
	  /* Set the pvt vector weights to 1 */
	  /* Not sure why this is needed since eahc particle has its own weight */

	  weights.push_back(1.0);
	  Particle temp_particle;
	  temp_particle.id = i;
	  temp_particle.x  = dist_x(gen);
	  temp_particle.y  = dist_y(gen);
	  temp_particle.theta= dist_theta(gen);
	  temp_particle.weight = 1.0;
	  particles.push_back(temp_particle);  
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/
        default_random_engine gen;	
	double sigma_x,sigma_y,sigma_theta;
	sigma_x = std_pos[0];
	sigma_y = std_pos[1];
	sigma_theta = std_pos[2];
	
	// Create distribution with Zero Mean
	// These are later added to final x/y/theta
//	normal_distribution<double> dist_x(0, sigma_x);
//	normal_distribution<double> dist_y(0, sigma_y);
//      normal_distribution<double> dist_theta(0, sigma_theta);

	for (int i=0;i<num_particles;i++) {
	  double current_x,current_y,current_theta;
	  double yr_dt;
	  double vel_by_yr;
	  current_x = particles[i].x;
	  current_y = particles[i].y;
	  current_theta = particles[i].theta;
	  yr_dt = yaw_rate * delta_t;
	  vel_by_yr = velocity/yaw_rate;

	  particles[i].x = current_x + (vel_by_yr * ( sin(current_theta + yr_dt)  - sin(current_theta)));
	  particles[i].y = current_y + (vel_by_yr * ( cos(current_theta ) - cos(current_theta + yr_dt)));
	  particles[i].theta = current_theta + yr_dt;

  	  normal_distribution<double> dist_x(particles[i].x, sigma_x);
    	  normal_distribution<double> dist_y(particles[i].y, sigma_y);
          normal_distribution<double> dist_theta(particles[i].theta, sigma_theta);

	  particles[i].x = dist_x(gen);
	  particles[i].y = dist_y(gen);
	  particles[i].theta = dist_theta(gen);

          #if 0
	  while (particles[i].theta > M_PI) {
	    particles[i].theta -= 2*M_PI;
	  }
	  while (particles[i].theta < -1*M_PI) {
	    particles[i].theta += 2*M_PI;
	  } //while
	  #endif
        
	} //i
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	// Steps  : For all observations, iterate through all map positions (predicted)
	// Find the minimum distance across all landmarks (map positions)
	// Update id of obersvation 

	for (int i=0;i<observations.size();i++) {
	  double minimum_distance = 10000.0;
	  for (int cnt =0;cnt<predicted.size();cnt++) {
	    double dist_calc; // Distance between the observation and map
	    dist_calc =  dist (observations[i].x,observations[i].y,predicted[cnt].x,predicted[cnt].y);
	    if (dist_calc < minimum_distance) {
	      /* Set new min distance and store id in observations */ 
	      minimum_distance = dist_calc; 
	      observations[i].id = predicted[cnt].id;
	    }
	  } // for cnt (predicted)
	
	} // for i
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	// std_landmark: Landmark measurment uncertainty
	// observations - landmark measurements
	// map_landmarks - map
	
        /* Before we start we clear the weights (global) vector */
	weights.clear();
        double sig_x,sig_y;
        sig_x = std_landmark[0];
        sig_y = std_landmark[1];

	for (int i=0; i<particles.size();i++) {
	// Step-1: For each map point take the particles positions and compute distance
	//         If the distance is within threshold, add to predicted vector
	
	  double x_part = particles[i].x;
	  double y_part = particles[i].y;
	  double theta  = particles[i].theta;

          double particle_weight = 1.0;
	  // Step1: For every observation, create the tranfomration to map
	  std::vector<LandmarkObs> observations_transf; /* Observations transformed to Map, based on particles position */
	  for (int cnt=0;cnt<observations.size();cnt++) { // Iterate through all observations and populate observations_transf
	    LandmarkObs obs_t; /* Temp for transforming every observation */
	    double x_obs,y_obs;

	    obs_t.id = observations[cnt].id;
	    x_obs    = observations[cnt].x;
	    y_obs    = observations[cnt].y;

	    obs_t.x  = x_part +  (cos(theta) * x_obs) - (sin(theta) * y_obs);
	    obs_t.y  = y_part +  (sin(theta) * x_obs) + (cos(theta) * y_obs); 
	    observations_transf.push_back(obs_t);

	    //cout << "Step-1 Completed for particle = "<<i<<endl;

            // Step-2: Make a landmark list depending on sensor range and distance from particle
	    // Iterate through the map, caculate distance of landmark from particle and if the
	    // distance is within sensor range, add the landmark to the predicted list
	    vector<LandmarkObs> predicted;
	    int id =0;
	    /* For debug, I ll use all landmark objects instead of filtered */
	    double minimum_distance = 1000.0;
	    double mu_x, mu_y;
            for (int map_id=0;map_id <map_landmarks.landmark_list.size();map_id++) {

              Map::single_landmark_s map_t; /* Temp for filtering landmarks */
              map_t = map_landmarks.landmark_list[map_id];
	      double dist_calc; // Distance between the observation and map
	      dist_calc = dist(obs_t.x,obs_t.y,map_t.x_f,map_t.y_f);
	      if (dist_calc < minimum_distance) { // New minima found
	        minimum_distance = dist_calc;
		mu_x =map_landmarks.landmark_list[map_id].x_f; 
		mu_y =map_landmarks.landmark_list[map_id].y_f; 
	      }
	    } // for map_id 

	    //cout << "Step-2 Completed for particle = "<<i<<endl;


	 
            // Step-3: Update associations, in this step each observation/transformed observation is
	    // mapped to a landmark
	    //cout << "Step-3 Completed for particle = "<<i<<endl;



	    // Step-4: Update weight

	    double weight_i; //Weight of the ith observation
	    x_obs  = obs_t.x;
	    y_obs  = obs_t.y;
	    double gauss_norm = 1.0/(2* M_PI * sig_x *sig_y);

	    double exponent   = (pow((x_obs - mu_x),2)/ (2 * sig_x * sig_x)) + (pow((y_obs - mu_y),2)/ (2 * sig_y * sig_y));

	    weight_i = gauss_norm * pow(2.71828,-1*exponent);

	    particle_weight = particle_weight * weight_i;

	    if (i==0) {
	      cout <<"theta" << theta<<"  x_obs = " <<observations[cnt].x << " "<< "mu_x ="<< mu_x <<"y_obs ="<<observations[cnt].y<<" mu_y = "<<mu_y<<endl;
	      cout <<"   "<< "gauss_norm ="<<gauss_norm<<"  "<< "exponent = "<<exponent << "  " << "weight_i = "<<weight_i<<endl;
	      cout <<"   "<< "particle_weight ="<<particle_weight<<endl;
	    }

          } // for cnt (observations)
	  particles[i].weight = particle_weight;

	  /* We also want to store the weights as a separate top level vector */
          weights.push_back(particle_weight);
	  if (i==0) {
	    cout << "Step-4 Completed for particle = "<<i<<endl;
	    cout << "Particle Weight = "<<particle_weight<<endl;
	  }

	} // for i
        #if 0
        cout <<"Weight of nth particle after update = " << particles[0].weight<<endl;
        cout <<"Weights of nth particle (global) after update = " << weights[0]<<endl;
        #endif
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> particles_new; // Create a new list ot populate

	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> d(weights.begin(), weights.end());
	std::map<int, int> m;
	for (int i=0;i<particles.size();i++) {
          particles_new.push_back(particles[d(gen)]); 
	}
        particles.clear(); /* Clean the old list and reassign */
	particles = particles_new;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
