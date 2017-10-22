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
from scipy.fftpack import fft

################################### Usr input part ##################################################

# 1 ##  point 3 #####################
x_target=0;
y_target=0;
z_upper_bound=14;
z_lower_bound=14;
tolerance=1e-3;
target_material_ID=1;
target_node_Index=0;


# 1 ##  point 3 #####################
# x_target=0;
# y_target=15;
# z_upper_bound=14;
# z_lower_bound=14;
# tolerance=1e-3;
# point_ID=3;
# target_material_ID=1;
# target_node_Index=0;


# # 2 ##  point 4 #####################
# x_target=0;
# y_target=15;
# z_upper_bound=0;
# z_lower_bound=0;
# tolerance=1e-3;
# point_ID=4;
# target_material_ID=1;
# target_node_Index=1;


# # 3 ##  point 5 #####################
# x_target=0;
# y_target=15;
# z_upper_bound=-36;
# z_lower_bound=-36;
# tolerance=1e-3;
# point_ID=5;
# target_material_ID=1;
# target_node_Index=2;


# # # 4 ##  point 1 #####################
# x_target=0;
# y_target=0;
# z_upper_bound=14;
# z_lower_bound=14;
# tolerance=1e-3;
# point_ID=1;
# target_material_ID=1;
# target_node_Index=3;


# ## 5 ## point 8 ###############
# x_target=0;
# y_target=15;
# z_upper_bound=0;
# z_lower_bound=0;
# tolerance=1e-3;
# point_ID=8;
# target_material_ID=2;
# target_node_Index=4;


# ## 6 ## point 9 ###############
# x_target=0;
# y_target=15;
# z_upper_bound=-36;
# z_lower_bound=-36;
# tolerance=1e-3;
# point_ID=9;
# target_material_ID=2;
# target_node_Index=5;


# ## 7 ##### point 12 ###############
# x_target=0;
# y_target=0;
# z_upper_bound=-36;
# z_lower_bound=-36;
# tolerance=1e-3;
# point_ID=12;
# target_material_ID=1;
# target_node_Index=6;


# ### 8 ##### point 13 ###############
# x_target=0;
# y_target=0;
# z_upper_bound=-36;
# z_lower_bound=-36;
# tolerance=1e-3;
# point_ID=13;
# target_material_ID=2;
# target_node_Index=7;





################# added this for depth variation trace plot [July,19,2017]###############################################
# x_target=0;

# y_target=15;

# tolerance=1e-3;



# target_material_ID=1;

target_node_Index=0;

##########################################################################################################################

No_mpi_process=4;


HDF5_file_prefix = "SMR_model2_SASSI_comparison_DRM_motion.h5"; 

# HDF5_file_prefix="SMR_DRM_analysis.h5";
# HDF5_file_prefix="free_field_model_DRM_analysis.h5";

################################# Ending usr input ###################################################


node=sp.loadtxt('data_analysis/node.txt')    ### node has five columns: one node-id; 3 node coordinates; one dof
element=sp.loadtxt('data_analysis/element.txt')  #### element has: one element-id; one element-type-tage(eg. 27 means 27-node-brick-element; 8 means 8-node-brick-element...); a bunch of node IDs that consisting of this element; one material ID;  

# print node.shape, element.shape
No_node=node.shape[0];
No_element=element.shape[0];
No_node_per_element=element.shape[1]-3;

target_node=sp.zeros((1,9));   ###Initialisze target node matrix has 9 columns:one target_node_index; one node ID, 3 node coordinates; 1 minimum element-ID that contains this node; materiaL-ID of element; one DOF; one displacement beginning row index in target node displacemnt matrix. 
target_node_component=sp.zeros((1,9));


target_node_displacement_index=1;


for x1 in xrange(0,No_node):
	######################## Depends on the geometrical description of traget nodes, this if condition need to be changed according to needs of user#######################
	# if abs(node[x1,1]-x_target)<tolerance and abs(node[x1,2]-y_target)<tolerance and abs(node[x1,3]-z_target)<tolerance:
	if abs(node[x1,1]-x_target)<tolerance and abs(node[x1,2]-y_target)<tolerance and node[x1,3]<=z_upper_bound and node[x1,3]>=z_lower_bound: 
	# if abs(node[x1,1]-x_target)<tolerance and abs(node[x1,2]-y_target)<tolerance and (abs(node[x1,3])<tolerance or abs(node[x1,3]+36)<tolerance or abs(node[x1,3]+8.394)<tolerance or abs(node[x1,3]+19.0147)<tolerance or abs(node[x1,3]+28.0275)<tolerance): 
	# if abs(node[x1,1]-x_target)<tolerance and abs(node[x1,2]-y_target)<tolerance and node[x1,3]<=0: 
	#######################################################################################################################################################################
		for x2 in xrange(0,No_element):
			for x3 in xrange(0,No_node_per_element):
				########################################## Depends on the geometrical description of traget nodes, this if condition need to be changed according to needs of user ####################################
				
				if element[x2,2+x3]==node[x1,0] and element[x2,No_node_per_element+2]==target_material_ID:

				# if element[x2,2+x3]==node[x1,0] and ((element[x2,No_node_per_element+2]==1) or (element[x2,No_node_per_element+2]==2)):
				
				# if element[x2,2+x3]==node[x1,0] and element[x2,No_node_per_element+2]==1:				
				#######################################################################################################################################################################################################
					target_node_Index=target_node_Index+1;
					target_node_component[0,0]=int(target_node_Index);                 ## target_node_ID (column 0)
					target_node_component[0,1]=int(node[x1,0]);                   ## corresponding ESSI node ID (column 1) 
					target_node_component[0,2:5]=node[x1,1:4]                     ## corresponding ESSI node coordinates (column 2,3,4).
					target_node_component[0,5]=int(element[x2,0]);                     ## element ID (column 5)
					target_node_component[0,6]=int(element[x2,No_node_per_element+2]); ## element material ID (column 6)
					target_node_component[0,7]=int(node[x1,4]);                        ## target_node_dof  (column 7)
					target_node_component[0,8]=int(target_node_displacement_index);    ## target_node_displacement beginning row index (column 8)
					target_node_displacement_index=target_node_displacement_index+node[x1,4];
					print "Target node No: ", int(target_node_component[0,0]), " at x=",target_node_component[0,2]," y=",target_node_component[0,3], " z=", target_node_component[0,4], " has been found at ESSI node:", int(target_node_component[0,1])," with ",int(target_node_component[0,7]),"dofs"," contained by element:", int(target_node_component[0,5]), "with material ID: ", int(target_node_component[0,6]),"\n";
					target_node=sp.concatenate((target_node,target_node_component));
					break;
			else:
				continue;
			break;
		# else:
		# 	continue;
		# break;

print "I am trying" ,target_node, target_node.shape

# master_file=HDF5_file_prefix+".feioutput";

master_file= "SMR_model2_SASSI_comparison_self_weight.h5.feioutput";   

f_master=h5py.File(master_file, "r");

node_partition=f_master['Model/Nodes/Partition'][:];   ## node_partition is column vector;

# Time=f_master['time'];    ## Time is a column vector;  Only one index is allowed. 

# time_step=Time[1]-Time[0];

# No_time_step=Time.shape[0];

#========Since DRM node partition info messed up, we get partition info from self weight file and get time vector from DRM motion file ===================
DRM_master_file= "SMR_model2_SASSI_comparison_DRM_motion.h5.feioutput";   

f_DRM_master=h5py.File(DRM_master_file, "r");

Time=f_DRM_master['time'];    ## Time is a column vector;  Only one index is allowed. 

time_step=Time[1]-Time[0];

No_time_step=Time.shape[0];

#========================================================================================================================================================

target_node_displacement=sp.zeros((1,No_time_step));              #### column number is timestep value. 

target_node_displacement_component=sp.zeros((1,No_time_step));

No_target_node=target_node.shape[0];

target_node_dof=0;

for x4 in xrange(0,No_target_node-1):	

	target_node_dof=int(target_node[x4+1,7]);

	target_node_partition=node_partition[int(target_node[x4+1,1])];

	print  "node partition is: ",  target_node_partition; 

	if target_node_partition<10:

		 target_hdf_file_postfix="."+"0"+str(target_node_partition)+".feioutput";  #target_hdf_file_postfix="."+str(target_node_partition)+".feioutput";  
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


################################################### Adding this part for calculation of acceleration, displacement_fft and acceleration_fft #####################################

total_target_node_DOFs=target_node_displacement.shape[0]-1;

target_node_acceleration=sp.zeros((total_target_node_DOFs,No_time_step-1));

target_node_acceleration_fft=sp.zeros((total_target_node_DOFs,No_time_step-1));

target_node_displacement_fft=sp.zeros((total_target_node_DOFs,No_time_step));

print "Number of time steps: ", No_time_step; 

print "Time steps: ",  time_step; 

l= range(0,No_time_step); 

# frequency=1.0/(No_time_step*time_step)*range(0,No_time_step);      ### sp.zeros((No_time_step,1));s
frequency=sp.array(l, dtype=int)*1.0/(No_time_step*time_step); 

print "I am printting frequency:\n ", frequency;   

for x6 in xrange(1,No_time_step-1):
		target_node_acceleration[:,x6]=(target_node_displacement[1:,x6+1]+target_node_displacement[1:,x6-1]-2*target_node_displacement[1:,x6])/(1.0*time_step*time_step);
		# print target_node_acceleration[:,x6];

for x7 in xrange(0,total_target_node_DOFs):
	target_node_displacement_fft[x7,:]=abs(fft(target_node_displacement[1+x7,:]))*2/No_time_step;
	target_node_acceleration_fft[x7,:]=abs(fft(target_node_acceleration[x7,:]))*2/(No_time_step-1);



# ######################################################### Added this part for find PGA along the depth #########################################################################
# ####################################################### Here we assume all target nodes are 3 DOFs######################################################################################

# No_rows_target_node_displacement=target_node_displacement.shape[0];

# # # print "I am verifying", target_node_displacement

# target_node_velocity=sp.zeros((No_rows_target_node_displacement,No_time_step));
# target_node_acceleration=sp.zeros((No_rows_target_node_displacement,No_time_step));

# for x6 in xrange(1,No_time_step):
# 	target_node_velocity[:,x6]=(target_node_displacement[:,x6]-target_node_displacement[:,x6-1])/time_step;

# for x7 in xrange(1,No_time_step):
# 	target_node_acceleration[:,x7]=(target_node_velocity[:,x7]-target_node_velocity[:,x7-1])/time_step;	


# abs_target_node_displacement=sp.absolute(target_node_displacement);

# abs_target_node_acceleration=sp.absolute(target_node_acceleration);

# pgu=sp.zeros((No_rows_target_node_displacement,1));  ##### This one is the composite one, contains displacements of all target nodes in three directions. ATTENTION: The first row is useless 

# pga=sp.zeros((No_rows_target_node_displacement,1));

# pgu=sp.amax(abs_target_node_displacement, axis=1);

# pga=sp.amax(abs_target_node_acceleration, axis=1);


# depth=sp.zeros((No_target_node,2));

# pgux=sp.zeros((No_target_node,1));
# pguy=sp.zeros((No_target_node,1));
# pguz=sp.zeros((No_target_node,1));

# pgax=sp.zeros((No_target_node,1));
# pgay=sp.zeros((No_target_node,1));
# pgaz=sp.zeros((No_target_node,1));

# for x8 in xrange(1,No_target_node):
# 	depth[x8,0]=target_node[x8,4];
# 	depth[x8,1]=x8;

# depth[0,0]=-200000000;  ###### Small trick to avoid the effects from first line

# depth=depth[sp.argsort(depth[:,0])]

# # print depth

# for x8 in xrange(1,No_target_node):
# 	index=depth[x8,1];
# 	pgux[x8]=pgu[1+3*(index-1)];
# 	pguy[x8]=pgu[2+3*(index-1)];
# 	pguz[x8]=pgu[3*index];
# 	pgax[x8]=pga[1+3*(index-1)];
# 	pgay[x8]=pga[2+3*(index-1)];
# 	pgaz[x8]=pga[3*index];



# plt.plot(pgux[1:,0],depth[1:,0])
# plt.xlabel('Peak displacement ' r'$U_x$ [m]')
# plt.ylabel('Depth Z [m]')
# plt.grid(True);
# output_fig_x="data_analysis/pgux.pdf";
# plt.axis('on')
# plt.savefig(output_fig_x)
# plt.show()

# plt.plot(pguy[1:,0],depth[1:,0])
# plt.xlabel('Peak displacement ' r'$U_y$ [m]')
# plt.ylabel('Depth Z [m]')
# plt.grid(True);
# output_fig_x="data_analysis/pguy.pdf";
# plt.axis('on')
# plt.savefig(output_fig_x)
# plt.show()

# plt.plot(pguz[1:,0],depth[1:,0])
# plt.xlabel('Peak displacement ' r'$U_z$ [m]')
# plt.ylabel('Depth Z [m]')
# plt.grid(True);
# output_fig_x="data_analysis/pguz.pdf";
# plt.axis('on')
# plt.savefig(output_fig_x)
# plt.show()


# plt.plot(0.1*pgax[1:,0],depth[1:,0])
# plt.xlabel('Peak acceleration ' r'$A_x$ [g]')
# plt.ylabel('Depth Z [m]')
# plt.grid(True);
# output_fig_x="data_analysis/pgax.pdf";
# plt.axis('on')
# plt.savefig(output_fig_x)
# plt.show()

# plt.plot(0.1*pgay[1:,0],depth[1:,0])
# plt.xlabel('Peak acceleration ' r'$A_y$ [g]')
# plt.ylabel('Depth Z [m]')
# plt.grid(True);
# output_fig_x="data_analysis/pgay.pdf";
# plt.axis('on')
# plt.savefig(output_fig_x)
# plt.show()

# plt.plot(0.1*pgaz[1:,0],depth[1:,0])
# plt.xlabel('Peak acceleration ' r'$A_z$ [g]')
# plt.ylabel('Depth Z [m]')
# plt.grid(True);
# output_fig_x="data_analysis/pgaz.pdf";
# plt.axis('on')
# plt.savefig(output_fig_x)
# plt.show()

# depth_filename='data_analysis/depth.txt';

# target_node_filename='data_analysis/target_node.txt';
# target_node_displacement_filename='data_analysis/target_node_displacement.txt';
# target_node_acceleration_filename='data_analysis/target_node_acceleration.txt';

# pgu_filename='data_analysis/pgu.txt';  # will have three columns: one column x; one column y and one columnz 

# pga_filename='data_analysis/pga.txt';  # will have three columns: one column x; one column y and one columnz 

# pgah_filename='data_analysis/pgah.txt';

# pguh_filename='data_analysis/pguh.txt';


# # pgux_filename='data_analysis/pgux.txt';
# # pguy_filename='data_analysis/pguy.txt';
# # pguz_filename='data_analysis/pguz.txt';

# # pgax_filename='data_analysis/pgax.txt';
# # pgay_filename='data_analysis/pgay.txt';
# # pgaz_filename='data_analysis/pgaz.txt';

# PGU=sp.zeros((No_target_node,1));
# PGA=sp.zeros((No_target_node,1));


# # PGU=[pgux, pguy, pguz];
# # PGA=[pgax, pgay, pgaz];


# PGU=sp.concatenate((PGU,pgux), axis=1);   ##### Please Note that effective data begins from rows number 1 and column number 1; row 0 and column 0 are not used. 
# PGU=sp.concatenate((PGU,pguy), axis=1);
# PGU=sp.concatenate((PGU,pguz), axis=1);

# PGA=sp.concatenate((PGA,pgax), axis=1);   ##### Please Note that effective data begins from rows number 1 and column number 1; row 0 and column 0 are not used. 
# PGA=sp.concatenate((PGA,pgay), axis=1);
# PGA=sp.concatenate((PGA,pgaz), axis=1);

# PGUH=sp.zeros((No_target_node,1));        ### Horizontal resultant displacement 
# PGAH=sp.zeros((No_target_node,1));		### Horizontal resultant acceleration

# for x9 in xrange(1,No_target_node):
# 	PGUH[x9]=pow(pow(PGU[x9][1],2)+pow(PGU[x9][2],2), 0.5);
# 	PGAH[x9]=pow(pow(PGA[x9][1],2)+pow(PGA[x9][2],2), 0.5);

# sp.savetxt(depth_filename,depth[1:,:]);

# sp.savetxt(target_node_filename,target_node);
# sp.savetxt(target_node_displacement_filename,target_node_displacement[1:,:]);
# sp.savetxt(target_node_acceleration_filename,target_node_acceleration[1:,:]);

# sp.savetxt(pgu_filename,PGU[1:,1:]);
# sp.savetxt(pga_filename,0.1*PGA[1:,1:]);


# sp.savetxt(pguh_filename,PGUH[1:]);
# sp.savetxt(pgah_filename,PGAH[1:]);

# # sp.savetxt(pgux_filename,pgux[1:]);
# # sp.savetxt(pguy_filename,pguy[1:]);
# # sp.savetxt(pguz_filename,pguz[1:]);

# # sp.savetxt(pgax_filename,0.1*pgax[1:]);
# # sp.savetxt(pgay_filename,0.1*pgay[1:]);
# # sp.savetxt(pgaz_filename,0.1*pgaz[1:]);

# sp.savetxt('data_analysis/Time.txt',Time);

# ##################################################### Ending find PGA part######################################################################################################






# ########################################## Usr can choose finish plotting here, or write another script, loading the output file and plot ######################################
plot_target_node_index=1;


# ########################################### This is filename applicable for extracting point by point ##############################################################################

target_node_filename='data_analysis/target_node_'+str(target_node[plot_target_node_index,0])+'.txt';
target_node_displacement_filename='data_analysis/target_node_displacement_'+str(target_node[plot_target_node_index,0])+'.txt';
target_node_acceleration_filename='data_analysis/target_node_acceleration_'+str(target_node[plot_target_node_index,0])+'.txt';
target_node_displacement_fft_filename='data_analysis/target_node_disfft_'+str(target_node[plot_target_node_index,0])+'.txt';
target_node_acceleration_fft_filename='data_analysis/target_node_accfft_'+str(target_node[plot_target_node_index,0])+'.txt';

# ###################################################################################################################################################################################

# ### OR

# ######################################### This is filename applicable for extracting many points at the same time ################################################################

# target_node_filename='data_analysis/target_node'+'.txt';
# target_node_displacement_filename='data_analysis/target_node_displacement'+'.txt';
# target_node_acceleration_filename='data_analysis/target_node_acceleration'+'.txt';
# target_node_displacement_fft_filename='data_analysis/target_node_disfft'+'.txt';
# target_node_acceleration_fft_filename='data_analysis/target_node_accfft'+'.txt';

# ###############################################################################################################################################################################






plot_target_node_dis_index=int(target_node[plot_target_node_index,8]);
ux=target_node_displacement[plot_target_node_dis_index,:];
uy=target_node_displacement[plot_target_node_dis_index+1,:];
uz=target_node_displacement[plot_target_node_dis_index+2,:];

ux=sp.transpose(ux);
uy=sp.transpose(uy);
uz=sp.transpose(uz);


plt.plot(Time,ux)
plt.xlabel('Time T / (s)')
plt.ylabel('Displacement U / (m)')
plt.title('Seismic Responce')
plt.grid(True);
# plt.axis([0.0,20.0,-0.25,0.05])
output_fig_x="data_analysis/point_"+str(target_node[plot_target_node_index,0])+"_x.pdf";
# plt.savefig(output_fig_x, transparent=True, bbox_inches='tight')
plt.axis('on')
plt.savefig(output_fig_x)
plt.show()



# # plt.subplot(3,1,2)
# plt.plot(Time,uy)
# plt.xlabel('Time T / (s)')
# plt.ylabel('Displacement V / (m)')
# plt.title('Seismic Responce')
# plt.grid()
# # plt.axis([0.0,20.0,-0.3,0.18])
# # plt.box()
# output_fig_y="data_analysis/point_"+str(target_node[plot_target_node_index,0])+"_y.pdf";
# plt.savefig(output_fig_y, transparent=True, bbox_inches='tight')
# plt.show()

# # plt.subplot(3,1,3)
# plt.plot(Time,uz)
# plt.xlabel('Time T / (s)')
# plt.ylabel('Displacement W / (m)')
# plt.title('Seismic Responce')
# plt.grid(True);
# # plt.axis([0.0,20.0,-0.12,0])
# # plt.box()
# output_fig_z="data_analysis/point_"+str(target_node[plot_target_node_index,0])+"_z.pdf";
# plt.savefig(output_fig_z, transparent=True, bbox_inches='tight')
# plt.show()

target_node_displacement=sp.transpose(target_node_displacement);
target_node_acceleration=sp.transpose(target_node_acceleration);
target_node_displacement_fft=sp.transpose(target_node_displacement_fft);
target_node_acceleration_fft=sp.transpose(target_node_acceleration_fft);
frequency=sp.transpose(frequency);

sp.savetxt(target_node_filename,target_node);
sp.savetxt(target_node_displacement_filename,target_node_displacement[:,1:]);
sp.savetxt('data_analysis/Time.txt',Time);
sp.savetxt('data_analysis/frequency.txt',frequency);
sp.savetxt(target_node_acceleration_filename,target_node_acceleration);
sp.savetxt(target_node_displacement_fft_filename,target_node_displacement_fft);
sp.savetxt(target_node_acceleration_fft_filename,target_node_acceleration_fft);

# ######################################################### Ending plot part #######################################################################