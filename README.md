# Global Peak Finding based on Distributed Gaussian Process and Active Sensing

This project is one part of the module: Intelligent Distributed System and was implemented by Weijie Qi and Yunru Qu under Prof. Dr. Daniele Fontanelli.

We implemented the solution for multiple robots working together to explore an unknown environment and eventually gather at the global peak. The process makes use of Distributed Gaussian Process and Average Consensus algorithms.

图（The GIF has been sped up for smoother display）

The system is designed to be distributed, where agents can only communicate with each other within a certain distance. Throughout the process, collision avoidance between agents and prevention of crossing map boundaries are taken into consideration.

![gif](https://github.com/Winnie-Qi/Global-Peak-Finding-based-on-Distributed-Gaussian-Process-and-Active-Sensing/mat/gif2.gif)

To run this project, please execute the main_mobile_agents.m file. You have the flexibility to modify the number of agents and the random seed to change the starting positions of the agents.

图（main_mobile_agents.m）

If you want to observe how information is updated and transferred between agents, please execute the main_stationary_sensors.m file. This script demonstrates the updating and transfer of information between stationary sensors. You can modify the density or the communication radius of the agents and choose which agent's results will be displayed.

图

图（main_stationary_sensors.m）

仓库具体的结构如下：