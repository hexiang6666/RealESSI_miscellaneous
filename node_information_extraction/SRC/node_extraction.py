#********************************************************************************************************
# File:              node_extraction.py                                              
# Author:            hexiang                                     | Boris Jeremic,                       
# Date:              2017-04-20 20:23:20                         | University of California, Davis,95616*
# Description:       #############                               | California                           #
# Rev:               Version 1                                   | jeremic@ucdavis.edu                  #
# Email:             hexwang@ucdavis.edu                         | Computational Geomechanics Group     #
# * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  # 
#                           Last Modified time: 2017-03-16 22:29:55                                     #              
#  * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #        
# The copyright to the computer program(s) herein is the property of Hexiang Wang and Boris Jeremic     #
# The program(s) may be used and/or copied only with written permission of Hexiang Wang or in accordance#
# with the terms and conditions stipulated in the agreement/contract under which the program have been  #
# supplied.                                                                                             #
#*******************************************************************************************************#

#! /usr/bin/env python
import scipy as sp
import h5py
import pickle
import matplotlib.pyplot as plt  

################################### Usr input part ##################################################
x_target=15;
y_target=15;
z_target=0;
tolerance=1e-3;
target_material_ID=2;

No_mpi_process=10;

HDF5_file_prefix="SMR_DRM_analysis.h5";

################################# Ending usr input ###################################################
node=sp.loadtxt('data_analysis/node.txt')    ### node has five columns: one node-id; 3 node coordinates; one dof
element=sp.loadtxt('data_analysis/element.txt')  #### element has: one element-id; one element-type-tage(eg. 27 means 27-node-brick-element; 8 means 8-node-brick-element...); a bunch of node IDs that consisting of this element; one material ID;  

# print node.shape, element.shape
No_node=node.shape[0];
No_element=element.shape[0];
No_node_per_element=element.shape[1]-3;

target_node=sp.zeros((1,9));   ###Initialisze target node matrix has 9 columns:one target_node_index; one node ID, 3 node coordinates; 1 minimum element-ID that contains this node; materiaL-ID of element; one DOF; one displacement beginning row index in target node displacemnt matrix. 
target_node_component=sp.zeros((1,9));

target_node_Index=0;
target_node_displacement_index=1;

for x1 in xrange(0,No_node):
	######################## Depends on the geometrical description of traget nodes, this if condition need to be changed according to needs of user#######################
	if abs(node[x1,1]-x_target)<tolerance and abs(node[x1,2]-y_target)<tolerance and abs(node[x1,3]-z_target)<tolerance:
	#######################################################################################################################################################################
		for x2 in xrange(0,No_element):
			for x3 in xrange(0,No_node_per_element):
				########################################## Depends on the geometrical description of traget nodes, this if condition need to be changed according to needs of user ####################################
				if element[x2,2+x3]==node[x1,0] and element[x2,No_node_per_element+2]==target_material_ID:
				#######################################################################################################################################################################################################
					target_node_Index=target_node_Index+1;
					target_node_component[0,0]=int(target_node_Index);                 ## target_node_ID (column 0)
					target_node_component[0,1]=int(node[x1,0]);                   ## corresponding ESSI node ID (column 1) 
					target_node_component[0,2:5]=node[x1,1:4]                     ## corresponding ESSI node coordinates (column 2,3,4).
					target_node_component[0,5]=int(element[x2,0]);                     ## element ID (column 5)
					target_node_component[0,6]=int(element[x2,No_node_per_element+2]); ## element material ID (column 6)
					target_node_component[0,7]=int(node[x1,4]);                        ## target_node_dof  (column 7)
					target_node_component[0,8]=int(target_node_displacement_index);    ## target_node_displacement beginning row index (column 8)
					target_node_displacement_Index=target_node_displacement_index+node[x1,4];
					print "Target node No: ", int(target_node_component[0,0]), " at x=",target_node_component[0,2]," y=",target_node_component[0,3], " z=", target_node_component[0,4], " has been found at ESSI node:", int(target_node_component[0,1])," with ",int(target_node_component[0,7]),"dofs"," contained by element:", int(target_node_component[0,5]), "with material ID: ", int(target_node_component[0,6]),"\n";
					target_node=sp.concatenate((target_node,target_node_component));
					break;
			else:
				continue;
			break;
		else:
			continue;
		break;

# print "I am trying" ,target_node 
master_file=HDF5_file_prefix+".feioutput";

f_master=h5py.File(master_file, "r");

node_partition=f_master['Model/Nodes/Partition'][:];   ## node_partition is column vector;

Time=f_master['time'];    ## Time is a column vector;

No_time_step=Time.shape[0];

# print node_partition.shape;

target_node_displacement=sp.zeros((1,No_time_step)); 
target_node_displacement_component=sp.zeros((1,No_time_step));

for x4 in xrange(0,target_node_Index):	
	target_node_dof=int(target_node[x4+1,7]);
	target_node_partition=node_partition[int(target_node[x4+1,1])];
	if target_node_partition<10:
		target_hdf_file_postfix="."+"0"+str(target_node_partition)+".feioutput";
	else:
		target_hdf_file_postfix="."+str(target_node_partition)+".feioutput";	
	target_hdf_file=HDF5_file_prefix+target_hdf_file_postfix;
	f_target=h5py.File(target_hdf_file,"r");
	Index_node_coordinates=f_target['Model/Nodes/Index_to_Coordinates'][int(target_node[x4+1,1])];
	node_coordinates=f_target['Model/Nodes/Coordinates'][:];
	Index_node_displacement=f_target['Model/Nodes/Index_to_Generalized_Displacements'][int(target_node[x4+1,1])];
	############# Double check if this node is out target node##################
	if abs(target_node[x4+1,2]-node_coordinates[Index_node_coordinates])>tolerance or abs(target_node[x4+1,3]-node_coordinates[Index_node_coordinates+1])>tolerance or abs(target_node[x4+1,4]-node_coordinates[Index_node_coordinates+2])>tolerance:
		print "Attention: find the wrong location in HDF5 for target node No:", x4+1, "Please check it!!!";
	else:
		for x5 in xrange(0,target_node_dof):
			target_node_displacement_component[0,:]=f_target['Model/Nodes/Generalized_Displacements'][Index_node_displacement+x5,:];
			target_node_displacement=sp.concatenate((target_node_displacement,target_node_displacement_component));   ##### Please Note that effective data begins from rows number 1; row 0 is not used


sp.savetxt('data_analysis/target_node.txt',target_node);
sp.savetxt('data_analysis/target_node_displacement.txt',target_node_displacement);
sp.savetxt('data_analysis/Time.txt',Time);


########################################### Usr can choose finish plotting here, or write another script, loading the output file and plot ######################################
plot_target_node_index=1;
plot_target_node_dis_index=int(target_node[plot_target_node_index,8]);
ux=target_node_displacement[plot_target_node_dis_index,:];
uy=target_node_displacement[plot_target_node_dis_index+1,:];
uz=target_node_displacement[plot_target_node_dis_index+2,:];

ux=sp.transpose(ux);
uy=sp.transpose(uy);
uz=sp.transpose(uz);


# fig_x=plt.figure();
# ax_x=fig_x.add_subplot(111)

# ax_x.axhline(y=-0.12,linewidth=1,color="black")  
# ax_x.axvline(linewidth=1,color="black")
# ax_x.axvline(x=20,linewidth=1,color="black")

plt.plot(Time,ux)
plt.xlabel('Time T / (s)')
plt.ylabel('Displacement U / (m)')
plt.title('Seismic Responce')
plt.grid(True);
plt.axis([0.0,20.0,-0.12,0.04])
output_fig_x="data_analysis/point_"+str(plot_target_node_index)+"_x.pdf";
# plt.savefig(output_fig_x, transparent=True, bbox_inches='tight')
plt.axis('on')
plt.savefig(output_fig_x)
plt.show()



# plt.subplot(3,1,2)
plt.plot(Time,uy)
plt.xlabel('Time T / (s)')
plt.ylabel('Displacement V / (m)')
plt.title('Seismic Responce')
plt.grid()
plt.axis([0.0,20.0,-0.02,0.1])
# plt.box()
output_fig_y="data_analysis/point_"+str(plot_target_node_index)+"_y.pdf";
plt.savefig(output_fig_y, transparent=True, bbox_inches='tight')
plt.show()

# plt.subplot(3,1,3)
plt.plot(Time,uz)
plt.xlabel('Time T / (s)')
plt.ylabel('Displacement W / (m)')
plt.title('Seismic Responce')
plt.grid(True);
plt.axis([0.0,20.0,-0.09,-0.01])
# plt.box()
output_fig_z="data_analysis/point_"+str(plot_target_node_index)+"_z.pdf";
plt.savefig(output_fig_z, transparent=True, bbox_inches='tight')
plt.show()