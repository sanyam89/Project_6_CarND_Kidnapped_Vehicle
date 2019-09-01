/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  //if (initialized)
  //  return;
  std::default_random_engine gen;
  num_particles = 25;  // TODO: Set the number of particles
  std::normal_distribution<> dist_x(x,std[0]);
  std::normal_distribution<> dist_y(y,std[1]);
  std::normal_distribution<> dist_theta(theta,std[2]);
  
  for (int i = 0; i<num_particles; ++i){
    Particle current_particle;
    current_particle.id = i;
    current_particle.x = dist_x(gen);
    current_particle.y = dist_y(gen);
    current_particle.theta = dist_theta(gen);
    current_particle.weight = 1.0;
    
    particles.push_back(current_particle);
  }
  //cout<<"Initialization complete "<<endl;
  //cout<<"real positions : x= "<<x <<" ; y = " <<y << " ; theta = "<<theta<<endl;
  //cout<<"dist positions : x= "<<particles[10].x <<" ; y = " <<particles[10].y << " ; theta = "<<particles[10].theta<<endl;
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  double x_new, y_new, theta_new;
  for (int i=0; i<num_particles; ++i){
    if (fabs(yaw_rate) < 0.001){
      x_new = particles[i].x + (velocity * delta_t * cos(particles[i].theta));
      y_new = particles[i].y + (velocity * delta_t * sin(particles[i].theta));
      theta_new = particles[i].theta;
    }
    else{
      theta_new = particles[i].theta + yaw_rate*delta_t;
      x_new = particles[i].x + (velocity/yaw_rate * (sin(theta_new) - sin(particles[i].theta)));
      y_new = particles[i].y + (velocity/yaw_rate * (cos(particles[i].theta) - cos(theta_new)));
    }
    
    std::normal_distribution<> dist_x(x_new,std_pos[0]);
    std::normal_distribution<> dist_y(y_new,std_pos[1]);
    std::normal_distribution<> dist_theta(theta_new,std_pos[2]);
    
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
  //cout<<"Prediction after vehicle movement complete "<<endl;
  //cout<<"predicted positions : x= "<<particles[10].x <<" ; y = " <<particles[10].y << " ; theta = "<<particles[10].theta<<endl;
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations, double sensor_range) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  unsigned int i, j;
  int closest_landmark_id;
  double current_distance, min_distance;
  for (i = 0; i<observations.size(); ++i){
    min_distance = sensor_range*2.0; // maximum value
    closest_landmark_id = -1;
    for(j = 0; j<predicted.size(); ++j){
      current_distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
      if (current_distance <= min_distance){
        min_distance = current_distance;
        closest_landmark_id = predicted[j].id;
      }
    }
    observations[i].id = closest_landmark_id;
    //cout<<"Closest landmark ID for "<< i<<" th observation = "<< closest_landmark_id << " with min dist = "<<min_distance<< endl;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double total_particles_weight = 0.0;
  for (int i=0; i<num_particles; ++i){
    // Transforming observations from vehicle coordinates to Map coordinates
    vector<LandmarkObs> transformed_observations;
    LandmarkObs transformed_obs;
    //cout<<"Particle "<<i<<"th position : x_part= "<<particles[i].x <<" ; y_part= " <<particles[i].y << " ; theta = "<<particles[i].theta<<endl;
    for (unsigned int j=0; j<observations.size(); ++j){
      transformed_obs.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
      transformed_obs.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
      transformed_obs.id = -1; // assign ID later using Data Association
      transformed_observations.push_back(transformed_obs);
      //cout<<"initial observation : x_obs = "<<observations[j].x <<" ; y_obs = " <<observations[j].y<<endl;
      //cout<<"Transformed observation : x_map = "<<transformed_obs.x <<" ; y_map = " <<transformed_obs.y<<endl;
  
    }
    // Find landmarks in range of the vehicle sensors
    vector<LandmarkObs> landmarks_in_range;
    double distance;
    Map::single_landmark_s current_landmark;
    //cout<<"total landmarks in map : " << map_landmarks.landmark_list.size()<<endl;
    
    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); ++j){
      current_landmark = map_landmarks.landmark_list[j];
      distance = dist(current_landmark.x_f,current_landmark.y_f, particles[i].x,particles[i].y);
      if (distance <= sensor_range*1.5)
        landmarks_in_range.push_back(LandmarkObs{int(j), current_landmark.x_f, current_landmark.y_f});
    }
    //cout<<"total landmarks in range of the vehicle : " << landmarks_in_range.size()<<endl;
    
    // Associate the nearest landmark with the transformed observations
    dataAssociation(landmarks_in_range,transformed_observations, sensor_range);
    //cout<<"Data association complete!"<<endl;
     double sigma_x = std_landmark[0];
     double sigma_y = std_landmark[1];
     double sigma_x2 = sigma_x*sigma_x;
     double sigma_y2 = sigma_y*sigma_y;
     double term1 = 1 / (2*M_PI*sigma_x*sigma_y);
     double term2;
    for (unsigned int j=0; j<transformed_observations.size(); ++j){
      for (unsigned int k=0; k<landmarks_in_range.size(); ++k){
        if (landmarks_in_range[k].id == transformed_observations[j].id){
          term2 = pow( transformed_observations[j].x - landmarks_in_range[k].x , 2 ) / (2*sigma_x2) + 
                  pow( transformed_observations[j].y - landmarks_in_range[k].y , 2 ) / (2*sigma_y2);
          particles[i].weight *= term1 * exp(-1*term2);
        }
      }
    }
  total_particles_weight += particles[i].weight;
  }
  for (int i=0; i<num_particles; ++i)
    particles[i].weight /=total_particles_weight;
}
 
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  double beta = 0.0;
  double max_weight = std::numeric_limits<double>::min();
  for (int i=0; i<num_particles; ++i)
    if (max_weight<particles[i].weight)
      max_weight = particles[i].weight;
  
  std::uniform_real_distribution <>dist_wt (0,max_weight*2.0);
  std::discrete_distribution <>rand_index (0,num_particles-1);
  vector<Particle> resampledParticles;
  int index = rand_index(gen);
  for (int i = 0; i< num_particles; ++i){
    beta += dist_wt(gen);
    while (beta > particles[index].weight){
      beta -= particles[index].weight;
      index = (index + 1 ) % num_particles;
    }
     resampledParticles.push_back(particles[index]);
  }
  particles = resampledParticles;  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}