#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Create symbols
        q, d, a, alpha = symbols("q d a alpha") # used for the generalized Modified DH Transformation matrix Tdh
        q1, q2, q3, q4, q5, q6, q7, q8 = symbols('q1:9')
        d1, d2, d3, d4, d5, d6, d7, d8 = symbols('d1:9')
        a0, a1, a2, a3, a4, a5, a6, a7 = symbols('a0:8')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7 = symbols('alpha0:8')
        
        # Create Modified DH parameters
        DHTable = { 	alpha0:	    0, a0:      0, d1:   0.75, 
        alpha1: -pi/2, a1:   0.35, d2:      0, q2: q2-pi/2,
        alpha2:     0, a2:   1.25, d3:      0,
        alpha3: -pi/2, a3: -0.054, d4:    1.5,
        alpha4:  pi/2, a4:      0, d5:      0,
        alpha5: -pi/2, a5:      0, d6:      0,
        alpha6:     0, a6:      0, d7: 0.303, q7: -pi/2,
        alpha7: -pi/2, a7:      0, d8:     0, q8: -pi/2 }	
        
        # Define Modified DH Transformation matrix
        Tdh = Matrix([ 	[ 	     cos(q), 	       -sin(q),		  0,		 a],
        [ sin(q)*cos(alpha), cos(q)*cos(alpha),	-sin(alpha), -sin(alpha)*d],
        [ sin(q)*sin(alpha), cos(q)*sin(alpha),	 cos(alpha),  cos(alpha)*d],
        [                 0,                 0,           0,             1]	])
        
        # Create individual transformation matrices
        T0_1 = Tdh.subs([(q, q1), (d, d1), (a,a0), (alpha, alpha0)])
        T0_1 = T0_1.subs(DHTable)
        T1_2 = Tdh.subs([(q, q2), (d, d2), (a,a1), (alpha, alpha1)])
        T1_2 = T1_2.subs(DHTable)
        T2_3 = Tdh.subs([(q, q3), (d, d3), (a,a2), (alpha, alpha2)])
        T2_3 = T2_3.subs(DHTable)
        T3_4 = Tdh.subs([(q, q4), (d, d4), (a,a3), (alpha, alpha3)])
        T3_4 = T3_4.subs(DHTable)
        T4_5 = Tdh.subs([(q, q5), (d, d5), (a,a4), (alpha, alpha4)])
        T4_5 = T4_5.subs(DHTable)
        T5_6 = Tdh.subs([(q, q6), (d, d6), (a,a5), (alpha, alpha5)])
        T5_6 = T5_6.subs(DHTable)
        T6_Gdh = Tdh.subs([(q, q7), (d, d7), (a,a6), (alpha, alpha6)])
        T6_Gdh = T6_Gdh.subs(DHTable)
        TGdh_Gurdf = Tdh.subs([(q, q8), (d, d8), (a,a7), (alpha, alpha7)])
        TGdh_Gurdf = TGdh_Gurdf.subs(DHTable)
        
        T0_2 = T0_1 * T1_2
        T0_3 = T0_2 * T2_3
        T0_4 = T0_3 * T3_4
        T0_5 = T0_4 * T4_5
        T0_6 = T0_5 * T5_6
        T0_Gdh = T0_6 * T6_Gdh
        T0_Gurdf = T0_Gdh * TGdh_Gurdf
        
        
        #### Inverse Kinematics
        
        ## steps to get the WC position w.r.t base frame
        # given EE pose => quaternion and EE Position w.r.t base frame
        # from quaternion get euler angles (ga, be, al)
        # from euler angles get rotation matix R0_Gurdftc
        # from R0_Gurdf and EE Position get Transformation matrix T0_Gurdftc
        # from T0_Gurdftc and T6_Gurdf get T0_6
        # from T0_6 get Position of WC w.r.t base frame
        al, be, ga = symbols("al be ga")
        R0_Gurdf = Matrix([	[  cos(al)*cos(be), cos(al)*sin(be)*sin(ga)-sin(al)*cos(ga), cos(al)*sin(be)*cos(ga)+sin(al)*sin(ga)],
        [           sin(al)*cos(be), sin(al)*sin(be)*sin(ga)+cos(al)*cos(ga),  sin(al)*sin(be)*cos(ga)-cos(al)*sin(ga)],
        [ -sin(be), cos(be)*sin(ga), cos(be)*cos(ga)] ])
        
        
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()
            
            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z
            
            # from quaternion get euler angles
            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
            [req.poses[x].orientation.x, req.poses[x].orientation.y,
            req.poses[x].orientation.z, req.poses[x].orientation.w])
            

            
            # from euler angles get rotation matix R0_Gurdftc
            R0_Gurdftc = R0_Gurdf.subs([(al, yaw), (be, pitch), (ga, roll) ])
            
            # from R0_Gurdftc and EE Position get Transformation matrix T0_Gurdftc
            T0_Gurdftc = R0_Gurdftc.col_insert(3, Matrix([px,py,pz]))
            T0_Gurdftc = T0_Gurdftc.row_insert(3, Matrix([[0,0,0,1]]))
            
            # get T6_Gurdf
            T6_Gurdf = T6_Gdh * TGdh_Gurdf
            
            # from T0_Gurdftc and T6_Gurdf get T0_6
            T0_6tc = T0_Gurdftc * T6_Gurdf.inv("LU")
            
            # get wc Position from T0_6
            Pwc = T0_6tc[:,3]
            
            # Calculate joint angles using Geometric IK method
            # the writeup illustrates how to get the joint angles geometrically
            theta1 = atan2(Pwc[1], Pwc[0])
            
            C = 1.25
            A = 1.50097168527591 # sqrt(0.054**2 + 1.5**2)
            D = Pwc[2]-0.75
            E = (sqrt(Pwc[0]**2 + Pwc[1]**2) - 0.35)
            B = sqrt(D**2 + E**2 )
            angA = acos( (A**2-B**2-C**2)/(-2*B*C) )
            angD = acos( (D**2-B**2-E**2)/(-2*B*E) )
            angB = acos( (B**2-A**2-C**2)/(-2*A*C) )
            theta2 = 1.570796327 - angA - angD
            
            theta3 = 1.570796327 - atan2(0.054, 1.5) - angB
            
            ## get equations from R3_6
            R0_3tc = T0_3[0:3, 0:3]
            R0_3tc = R0_3tc.evalf(subs={q1: theta1, q2:theta2, q3:theta3})
            R3_6tc = R0_3tc.inv("LU") * T0_6tc[0:3,0:3]
            theta4 = atan2(R3_6tc[2,2], -R3_6tc[0,2])
            theta5 = atan2(sqrt( R3_6tc[0,2]**2 + R3_6tc[2,2]**2 ), R3_6tc[1,2])
            theta6 = atan2(-R3_6tc[1,1], R3_6tc[1,0])
            ###
            
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)
        
        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
