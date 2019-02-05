#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{
	auto found_null = s.find("null");
	auto b1 = s.find_first_of("[");
	auto b2 = s.find_first_of("}");
	if (found_null != string::npos)
	{
		return "";
	}
	else if (b1 != string::npos && b2 != string::npos)
	{
		return s.substr(b1, b2 - b1 + 2);
	}
	return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for (int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x, y, map_x, map_y);
		if (dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}

	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y - y), (map_x - x));

	double angle = fabs(theta - heading);
	angle = min(2 * pi() - angle, angle);

	if (angle > pi() / 2)
	{
		closestWaypoint++;
		if (closestWaypoint == maps_x.size())
		{
			closestWaypoint = 0;
		}
	}

	return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp - 1;
	if (next_wp == 0)
	{
		prev_wp = maps_x.size() - 1;
	}

	double n_x = maps_x[next_wp] - maps_x[prev_wp];
	double n_y = maps_y[next_wp] - maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
	double proj_x = proj_norm * n_x;
	double proj_y = proj_norm * n_y;

	double frenet_d = distance(x_x, x_y, proj_x, proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x, center_y, x_x, x_y);
	double centerToRef = distance(center_x, center_y, proj_x, proj_y);

	if (centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for (int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
	}

	frenet_s += distance(0, 0, proj_x, proj_y);

	return {frenet_s, frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp + 1) % maps_x.size();

	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
	double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

	double perp_heading = heading - pi() / 2;

	double x = seg_x + d * cos(perp_heading);
	double y = seg_y + d * sin(perp_heading);

	return {x, y};
}

int main()
{
	uWS::Hub h;

	// Load up map values for waypoint's x,y,s and d normalized normal vectors
	vector<double> map_waypoints_x;
	vector<double> map_waypoints_y;
	vector<double> map_waypoints_s;
	vector<double> map_waypoints_dx;
	vector<double> map_waypoints_dy;

	// Waypoint map to read from
	string map_file_ = "../data/highway_map.csv";
	// The max s value before wrapping around the track back to 0
	double max_s = 6945.554;

	ifstream in_map_(map_file_.c_str(), ifstream::in);

	string line;
	while (getline(in_map_, line))
	{
		istringstream iss(line);
		double x;
		double y;
		float s;
		float d_x;
		float d_y;
		iss >> x;
		iss >> y;
		iss >> s;
		iss >> d_x;
		iss >> d_y;
		map_waypoints_x.push_back(x);
		map_waypoints_y.push_back(y);
		map_waypoints_s.push_back(s);
		map_waypoints_dx.push_back(d_x);
		map_waypoints_dy.push_back(d_y);
	}

	h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
																											 uWS::OpCode opCode) {
		// "42" at the start of the message means there's a websocket message event.
		// The 4 signifies a websocket message
		// The 2 signifies a websocket event
		//auto sdata = string(data).substr(0, length);
		//cout << sdata << endl;
		if (length && length > 2 && data[0] == '4' && data[1] == '2')
		{

			auto s = hasData(data);

			if (s != "")
			{
				auto j = json::parse(s);

				string event = j[0].get<string>();

				if (event == "telemetry")
				{
					// j[1] is the data JSON object

					// Main car's localization Data
					double car_x = j[1]["x"];
					double car_y = j[1]["y"];
					double car_s = j[1]["s"];
					double car_d = j[1]["d"];
					double car_yaw = j[1]["yaw"];
					double car_speed = j[1]["speed"];

					// Previous path data given to the Planner
					auto previous_path_x = j[1]["previous_path_x"];
					auto previous_path_y = j[1]["previous_path_y"];
					// Previous path's end s and d values
					double end_path_s = j[1]["end_path_s"];
					double end_path_d = j[1]["end_path_d"];

					// Sensor Fusion Data, a list of all other cars on the same side of the road.
					auto sensor_fusion = j[1]["sensor_fusion"];

					json msgJson;

					double target_speed = 22; //m/s
					double MAX_ACC = 5.0;	  //m/s^2
					double INTERVAL = 20.0 /1000.0; // 20ms
					double BUFFER_DISTANCE = 22; //m

					vector<double> next_x_vals;
					vector<double> next_y_vals;


					//cout << "Next Speed: " <<  next_speed << ", Distance POints:" << distance_betwee_points << endl;

					// for(int i =0 ; i< 50; i++){
					// 	next_x_vals.push_back(car_x + distance_betwee_points* i * cos(deg2rad(car_yaw)));
					// 	next_y_vals.push_back(car_y + distance_betwee_points* i * sin(deg2rad(car_yaw)));
					// }
					double pos_x;
					double pos_y;
					double angle;
					double last_speed = 0;
					int path_size = previous_path_x.size();
					for (int i = 0; i < path_size; ++i)
					{
						next_x_vals.push_back(previous_path_x[i]);
						next_y_vals.push_back(previous_path_y[i]);
					}
					if (path_size == 0)
					{
						pos_x = car_x;
						pos_y = car_y;
						angle = deg2rad(car_yaw);
						last_speed = car_speed;
					}
					else
					{
						pos_x = previous_path_x[path_size - 1];
						pos_y = previous_path_y[path_size - 1];

						
						
						double pos_x2 = previous_path_x[path_size - 2];
						double pos_y2 = previous_path_y[path_size - 2];
						angle = atan2(pos_y - pos_y2, pos_x - pos_x2);
						last_speed = sqrt ((pos_x2 - pos_x) * (pos_x2 - pos_x) + (pos_y2-pos_y) * (pos_y2-pos_y))/0.02;
					}

					//Add first point of previous path back to waypoints for spline
					vector<double> s_waypoints_x;
					vector<double> s_waypoints_y;
					
					s_waypoints_x.push_back(pos_x);
					s_waypoints_y.push_back(pos_y);


					//Drive along the road
					//Find s & d of last known points
					//cout << end_path_s << ", " << end_path_d << endl;
					cout << "POS: " << pos_x << ", " << pos_y << "-" << path_size << "," << angle << endl;
					vector<double> sd = getFrenet(pos_x, pos_y, angle, map_waypoints_x, map_waypoints_y);
						cout << sd[0] << " ," << sd[1] << endl;
					
					
					// Localize other cars
					bool car_front = false;
					bool car_left = false;
					bool car_right = false;

					//What is my current lane?
					int ego_lane = -1;
					if(sd[1] > 0 && sd[1] <=4){
						ego_lane = 0;
					}else if(sd[1] >4 && sd[1] <=8){
						ego_lane = 1;
					}else if(sd[1]>8 && sd[1] <=12){
						ego_lane =2;
					}else{
						cout << "Ego out of mind (d): " << sd[1] << endl;
						exit(0);
					}
					//cout << "Ego Lane: " << ego_lane << endl;

					//Is anyone infront (assume 2 seconds look ahead. 22m/s *2s = 44m)
					//cout << "Num cars: "  << sensor_fusion.size() << endl;
					for(int i=0;i< sensor_fusion.size(); i++){
						double car_d = sensor_fusion[i][6];
						int car_lane = -1;
						if(car_d >0 && car_d <=4){
							car_lane = 0;	
						}else if(car_d > 4 && car_d <=8){
							car_lane = 1;
						}else if(car_d >8 && car_d <=12){
							car_lane = 2;
						}
						if(car_lane == -1){
							cout << "Unknown lane car: " << sensor_fusion[i] << endl;
							continue;
						}
						if(car_lane == ego_lane){
							double s_ego = sd[0];
							double s_car = sensor_fusion[i][5];
							//cout << "Car Infront at distance" << s_car - s_ego<< endl;
							double distance_front = (s_car - s_ego);
							if( distance_front < BUFFER_DISTANCE && distance_front > 0){
								
								//cout << sensor_fusion[i] << endl;
								// Follow car 
								double vx = sensor_fusion[i][3];
								double vy = sensor_fusion[i][4];
								double car_v = sqrt(vx*vx + vy*vy);
								cout << "Car in range (m):" << ( s_car - s_ego) << " ,v:"<< car_v <<  endl;
								target_speed = car_v;
							}
						}
					}

					double s_increment = sd[0];

					//Get rough points for spline connection (for going straight)
					double max_s = 100; //2 * (car_speed * 1 + 0.5 * MAX_ACC * 1 * 1) ; //max distance at MAX_ACC x 2
					int NUM_SPLINE_POINTS = 5;
					double s_spline_increment = (max_s)/double(NUM_SPLINE_POINTS);
					
					for (int i=1;i< NUM_SPLINE_POINTS; i++){
						vector<double> xy = getXY( s_increment + i*s_spline_increment, 6.0, map_waypoints_s, map_waypoints_x, map_waypoints_y);	
						//cout << "Spline: " << xy[0] << ", "<< xy[1] << endl;
						s_waypoints_x.push_back(xy[0]);
						s_waypoints_y.push_back(xy[1]);
					}

					// If we use these x,y coordinates directly, then there is no guarantee that x is linearly increasing
					// since the map is in real world coordinates. Hence we need to convert to car coordinates

					//Convert from real world to car
					for (int i=0; i< s_waypoints_x.size(); i++ ){
						double shift_x = s_waypoints_x[i] - pos_x;
						double shift_y = s_waypoints_y[i] - pos_y;

						s_waypoints_x[i] = shift_x * cos (-angle) - shift_y * sin(-angle);
						s_waypoints_y[i] = shift_x * sin (-angle) + shift_y * cos(-angle);
					}

					//Apply the spline
					tk::spline s;
   					s.set_points(s_waypoints_x,s_waypoints_y);
					
					// Drive with linear acceleration along the spline starting at zero
					double start_origin_x = 0;
					for (int i = 0; i < 50 - path_size; ++i)
					{					
						
						double ACC = MAX_ACC;
						// If target speed is not reached, keep accelerating else acceleration = 0
						if (last_speed > target_speed)
						{
							last_speed = last_speed - MAX_ACC * INTERVAL;
							ACC = - 2 * MAX_ACC;
						} else if (last_speed == target_speed){
							last_speed = target_speed;
							ACC = 0.0;
						} else if (last_speed < target_speed){
							last_speed = last_speed + MAX_ACC * INTERVAL;
							ACC = MAX_ACC;
						}
						double distance_betwee_points = last_speed * INTERVAL + 0.5 * ACC * INTERVAL * INTERVAL;
						//cout << next_speed << ", " << distance_betwee_points << endl;
						double x_point = start_origin_x + distance_betwee_points;
						double y_point = s(x_point);

						start_origin_x = x_point;
						
						double x_ref = x_point;
						double y_ref = y_point;

						x_point = x_ref * cos(angle) - y_ref * sin(angle);
              			y_point = x_ref * sin(angle) + y_ref * cos(angle);

						x_point += pos_x;
						y_point += pos_y;

						//cout << x_point << "," << y_point<< endl;
						

						next_x_vals.push_back(x_point);
						next_y_vals.push_back(y_point);
					}
					//exit(0);
					// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
					msgJson["next_x"] = next_x_vals;
					msgJson["next_y"] = next_y_vals;

					auto msg = "42[\"control\"," + msgJson.dump() + "]";

					//this_thread::sleep_for(chrono::milliseconds(1000));
					ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
				}
			}
			else
			{
				// Manual driving
				std::string msg = "42[\"manual\",{}]";
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			}
		}
	});

	// We don't need this since we're not using HTTP but if it's removed the
	// program
	// doesn't compile :-(
	h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
					   size_t, size_t) {
		const std::string s = "<h1>Hello world!</h1>";
		if (req.getUrl().valueLength == 1)
		{
			res->end(s.data(), s.length());
		}
		else
		{
			// i guess this should be done more gracefully?
			res->end(nullptr, 0);
		}
	});

	h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
		std::cout << "Connected!!!" << std::endl;
	});

	h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
						   char *message, size_t length) {
		ws.close();
		std::cout << "Disconnected" << std::endl;
	});

	int port = 4567;
	if (h.listen(port))
	{
		std::cout << "Listening to port " << port << std::endl;
	}
	else
	{
		std::cerr << "Failed to listen to port" << std::endl;
		return -1;
	}
	h.run();
}
