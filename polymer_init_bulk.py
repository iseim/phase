def polymer_init_bulk(RNA_num, RNA_len, pro_num, pro_len, extent):

	import numpy as np
	from scipy.spatial import distance

	#initialize output position list
	poslist = list();

	#initialize position lists for each octant (only used internally)
	poslist1 = list(); #x>0, y>0, z>0
	poslist2 = list(); #x<0, y>0, z>0
	poslist3 = list(); #x<0, y<0, z>0
	poslist4 = list(); #x>0, y<0, z>0
	poslist5 = list(); #x>0, y>0, z<0
	poslist6 = list(); #x<0, y>0, z<0
	poslist7 = list(); #x<0, y<0, z<0
	poslist8 = list(); #x>0, y<0, z<0


	#initial position of first particle
	initx = np.random.uniform(low=-extent/2, high=extent/2);
	inity = np.random.uniform(low=-extent/2, high=extent/2);
	initz = np.random.uniform(low=-extent/2, high=extent/2);

	poslist.append([initx, inity, initz]);

	if initx>0 and inity>0 and initz>0:
		poslist1.append([initx, inity, initz]);
	elif initx<0 and inity>0 and initz>0:
		poslist2.append([initx, inity, initz]);
	elif initx<0 and inity<0 and initz>0:
		poslist3.append([initx, inity, initz]);
	elif initx>0 and inity<0 and initz>0:
		poslist4.append([initx, inity, initz]);
	elif initx>0 and inity>0 and initz<0:
		poslist5.append([initx, inity, initz]);
	elif initx<0 and inity>0 and initz<0:
		poslist6.append([initx, inity, initz]);
	elif initx<0 and inity<0 and initz<0:
		poslist7.append([initx, inity, initz]);
	elif initx>0 and inity<0 and initz<0:
		poslist8.append([initx, inity, initz]);


	index = 0;

	#loop through RNAs first
	for r in range(RNA_num):

		if len(poslist) > 1:

			while True:

				#candidate initial position of new RNA
				initx = np.random.uniform(low=-extent/2, high=extent/2);
				inity = np.random.uniform(low=-extent/2, high=extent/2);
				initz = np.random.uniform(low=-extent/2, high=extent/2);
				point = np.array([[initx, inity, initz]]);


				#check which octant candidate particle occupies and check for overlaps in corresponding list

				#octant 1
				if initx>0 and inity>0 and initz>0:
					if len(poslist1) == 0:
						check2 = False; check4 = False; check5 = False;
						if initx < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check2) or np.any(check4) or np.any(check5)):
							poslist.append([initx, inity, initz]);
							poslist1.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist1_array = np.asarray(poslist1);
						dist = distance.cdist(poslist1_array, point);
						check = np.less(dist, 1); check2 = False; check4 = False; check5 = False;
						if initx < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check2) or np.any(check4) or np.any(check5)):
							poslist.append([initx, inity, initz]);
							poslist1.append([initx, inity, initz]);
							index += 1
							break

				#octant 2
				elif initx<0 and inity>0 and initz>0:
					if len(poslist2) == 0:
						check1 = False; check3 = False; check6 = False;
						if initx > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check1) or np.any(check3) or np.any(check6)):
							poslist.append([initx, inity, initz]);
							poslist2.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist2_array = np.asarray(poslist2);
						dist = distance.cdist(poslist2_array, point);
						check = np.less(dist, 1); check1 = False; check3 = False; check6 = False;
						if initx > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check1) or np.any(check3) or np.any(check6)):
							poslist.append([initx, inity, initz]);
							poslist2.append([initx, inity, initz]);
							index += 1
							break

				#octant 3
				elif initx<0 and inity<0 and initz>0:
					if len(poslist3) == 0:
						check4 = False; check2 = False; check7 = False;
						if initx > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check4) or np.any(check2) or np.any(check7)):
							poslist.append([initx, inity, initz]);
							poslist3.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist3_array = np.asarray(poslist3);
						dist = distance.cdist(poslist3_array, point);
						check = np.less(dist, 1); check4 = False; check2 = False; check7 = False;
						if initx > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check4) or np.any(check2) or np.any(check7)):
							poslist.append([initx, inity, initz]);
							poslist3.append([initx, inity, initz]);
							index += 1
							break

				#octant 4
				elif initx>0 and inity<0 and initz>0:
					if len(poslist4) == 0:
						check3 = False; check1 = False; check8 = False;
						if initx < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check3) or np.any(check1) or np.any(check8)):
							poslist.append([initx, inity, initz]);
							poslist4.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist4_array = np.asarray(poslist4);
						dist = distance.cdist(poslist4_array, point);
						check = np.less(dist, 1); check3 = False; check1 = False; check8 = False;
						if initx < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check3) or np.any(check1) or np.any(check8)):
							poslist.append([initx, inity, initz]);
							poslist4.append([initx, inity, initz]);
							index += 1
							break

				#octant 5
				elif initx>0 and inity>0 and initz<0:
					if len(poslist5) == 0:
						check6 = False; check8 = False; check1 = False;
						if initx < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check6) or np.any(check8) or np.any(check1)):
							poslist.append([initx, inity, initz]);
							poslist5.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist5_array = np.asarray(poslist5);
						dist = distance.cdist(poslist5_array, point);
						check = np.less(dist, 1); check6 = False; check8 = False; check1 = False;
						if initx < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check6) or np.any(check8) or np.any(check1)):
							poslist.append([initx, inity, initz]);
							poslist5.append([initx, inity, initz]);
							index += 1
							break

				#octant 6
				elif initx<0 and inity>0 and initz<0:
					if len(poslist6) == 0:
						check5 = False; check7 = False; check2 = False;
						if initx > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check5) or np.any(check7) or np.any(check2)):
							poslist.append([initx, inity, initz]);
							poslist6.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist6_array = np.asarray(poslist6);
						dist = distance.cdist(poslist6_array, point);
						check = np.less(dist, 1); check5 = False; check7 = False; check2 = False;
						if initx > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check5) or np.any(check7) or np.any(check2)):
							poslist.append([initx, inity, initz]);
							poslist6.append([initx, inity, initz]);
							index += 1
							break

				#octant 7
				elif initx<0 and inity<0 and initz<0:
					if len(poslist7) == 0:
						check8 = False; check6 = False; check3 = False;
						if initx > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check8) or np.any(check6) or np.any(check3)):
							poslist.append([initx, inity, initz]);
							poslist7.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist7_array = np.asarray(poslist7);
						dist = distance.cdist(poslist7_array, point);
						check = np.less(dist, 1); check8 = False; check6 = False; check3 = False;
						if initx > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check8) or np.any(check6) or np.any(check3)):
							poslist.append([initx, inity, initz]);
							poslist7.append([initx, inity, initz]);
							index += 1
							break

				#octant 8
				elif initx>0 and inity<0 and initz<0:
					if len(poslist8) == 0:
						check7 = False; check5 = False; check4 = False;
						if initx < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check7) or np.any(check5) or np.any(check4)):
							poslist.append([initx, inity, initz]);
							poslist8.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist8_array = np.asarray(poslist8);
						dist = distance.cdist(poslist8_array, point);
						check = np.less(dist, 1); check7 = False; check5 = False; check4 = False;
						if initx < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check7) or np.any(check5) or np.any(check4)):
							poslist.append([initx, inity, initz]);
							poslist8.append([initx, inity, initz]);
							index += 1
							break


		#build up the rest of the RNA
		for i in range(RNA_len-1):

			while True:

				#candidate position of new particle
				pi = np.random.rand(1) * 2*np.pi;
				z_it = np.random.uniform(low=-1, high=1);
				x = np.sqrt(1-(z_it**2))*np.cos(pi) + poslist[index][0];
				y = np.sqrt(1-(z_it**2))*np.sin(pi) + poslist[index][1];
				z = z_it + poslist[index][2];
				point = np.array([[x, y, z]]);


				#check for overlaps

				#octant 1
				if x>0 and y>0 and z>0:
					if len(poslist1) == 0:
						check2 = False; check4 = False; check5 = False;
						if x < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check2) or np.any(check4) or np.any(check5)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist1.append([x, y, z]);
								index += 1
								break
					else:
						poslist1_array = np.asarray(poslist1);
						dist = distance.cdist(poslist1_array, point);
						check = np.less(dist, 1); check2 = False; check4 = False; check5 = False;
						if x < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check2) or np.any(check4) or np.any(check5)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist1.append([x, y, z]);
								index += 1
								break

				#octant 2
				elif x<0 and y>0 and z>0:
					if len(poslist2) == 0:
						check1 = False; check3 = False; check6 = False;
						if x > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check1) or np.any(check3) or np.any(check6)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist2.append([x, y, z]);
								index += 1
								break
					else:
						poslist2_array = np.asarray(poslist2);
						dist = distance.cdist(poslist2_array, point);
						check = np.less(dist, 1); check1 = False; check3 = False; check6 = False;
						if x > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check1) or np.any(check3) or np.any(check6)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist2.append([x, y, z]);
								index += 1
								break

				#octant 3
				elif x<0 and y<0 and z>0:
					if len(poslist3) == 0:
						check4 = False; check2 = False; check7 = False;
						if x > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check4) or np.any(check2) or np.any(check7)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist3.append([x, y, z]);
								index += 1
								break
					else:
						poslist3_array = np.asarray(poslist3);
						dist = distance.cdist(poslist3_array, point);
						check = np.less(dist, 1); check4 = False; check2 = False; check7 = False;
						if x > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check4) or np.any(check2) or np.any(check7)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist3.append([x, y, z]);
								index += 1
								break

				#octant 4
				elif x>0 and y<0 and z>0:
					if len(poslist4) == 0:
						check3 = False; check1 = False; check8 = False;
						if x < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check3) or np.any(check1) or np.any(check8)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist4.append([x, y, z]);
								index += 1
								break
					else:
						poslist4_array = np.asarray(poslist4);
						dist = distance.cdist(poslist4_array, point);
						check = np.less(dist, 1); check3 = False; check1 = False; check8 = False;
						if x < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check3) or np.any(check1) or np.any(check8)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist4.append([x, y, z]);
								index += 1
								break

				#octant 5
				elif x>0 and y>0 and z<0:
					if len(poslist5) == 0:
						check6 = False; check8 = False; check1 = False;
						if x < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check6) or np.any(check8) or np.any(check1)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist5.append([x, y, z]);
								index += 1
								break
					else:
						poslist5_array = np.asarray(poslist5);
						dist = distance.cdist(poslist5_array, point);
						check = np.less(dist, 1); check6 = False; check8 = False; check1 = False;
						if x < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check6) or np.any(check8) or np.any(check1)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist5.append([x, y, z]);
								index += 1
								break

				#octant 6
				elif x<0 and y>0 and z<0:
					if len(poslist6) == 0:
						check5 = False; check7 = False; check2 = False;
						if x > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check5) or np.any(check7) or np.any(check2)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist6.append([x, y, z]);
								index += 1
								break
					else:
						poslist6_array = np.asarray(poslist6);
						dist = distance.cdist(poslist6_array, point);
						check = np.less(dist, 1); check5 = False; check7 = False; check2 = False;
						if x > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check5) or np.any(check7) or np.any(check2)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist6.append([x, y, z]);
								index += 1
								break

				#octant 7
				elif x<0 and y<0 and z<0:
					if len(poslist7) == 0:
						check8 = False; check6 = False; check3 = False;
						if x > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check8) or np.any(check6) or np.any(check3)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist7.append([x, y, z]);
								index += 1
								break
					else:
						poslist7_array = np.asarray(poslist7);
						dist = distance.cdist(poslist7_array, point);
						check = np.less(dist, 1); check8 = False; check6 = False; check3 = False;
						if x > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check8) or np.any(check6) or np.any(check3)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist7.append([x, y, z]);
								index += 1
								break

				#octant 8
				elif x>0 and y<0 and z<0:
					if len(poslist8) == 0:
						check7 = False; check5 = False; check4 = False;
						if x < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check7) or np.any(check5) or np.any(check4)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist8.append([x, y, z]);
								index += 1
								break
					else:
						poslist8_array = np.asarray(poslist8);
						dist = distance.cdist(poslist8_array, point);
						check = np.less(dist, 1); check7 = False; check5 = False; check4 = False;
						if x < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check7) or np.any(check5) or np.any(check4)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist8.append([x, y, z]);
								index += 1
								break



	#loop through proteins
	for r in range(pro_num):

		if len(poslist) > 1 and RNA_num > 0 or len(poslist) == 1:

			while True:

				#candidate initial position of new protein
				initx = np.random.uniform(low=-extent/2, high=extent/2);
				inity = np.random.uniform(low=-extent/2, high=extent/2);
				initz = np.random.uniform(low=-extent/2, high=extent/2);
				point = np.array([[initx, inity, initz]]);

				#check which octant candidate particle occupies and check for overlaps in corresponding list
				#octant 1
				if initx>0 and inity>0 and initz>0:
					if len(poslist1) == 0:
						check2 = False; check4 = False; check5 = False;
						if initx < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check2) or np.any(check4) or np.any(check5)):
							poslist.append([initx, inity, initz]);
							poslist1.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist1_array = np.asarray(poslist1);
						dist = distance.cdist(poslist1_array, point);
						check = np.less(dist, 1); check2 = False; check4 = False; check5 = False;
						if initx < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check2) or np.any(check4) or np.any(check5)):
							poslist.append([initx, inity, initz]);
							poslist1.append([initx, inity, initz]);
							index += 1
							break

				#octant 2
				elif initx<0 and inity>0 and initz>0:
					if len(poslist2) == 0:
						check1 = False; check3 = False; check6 = False;
						if initx > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check1) or np.any(check3) or np.any(check6)):
							poslist.append([initx, inity, initz]);
							poslist2.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist2_array = np.asarray(poslist2);
						dist = distance.cdist(poslist2_array, point);
						check = np.less(dist, 1); check1 = False; check3 = False; check6 = False;
						if initx > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check1) or np.any(check3) or np.any(check6)):
							poslist.append([initx, inity, initz]);
							poslist2.append([initx, inity, initz]);
							index += 1
							break

				#octant 3
				elif initx<0 and inity<0 and initz>0:
					if len(poslist3) == 0:
						check4 = False; check2 = False; check7 = False;
						if initx > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check4) or np.any(check2) or np.any(check7)):
							poslist.append([initx, inity, initz]);
							poslist3.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist3_array = np.asarray(poslist3);
						dist = distance.cdist(poslist3_array, point);
						check = np.less(dist, 1); check4 = False; check2 = False; check7 = False;
						if initx > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check4) or np.any(check2) or np.any(check7)):
							poslist.append([initx, inity, initz]);
							poslist3.append([initx, inity, initz]);
							index += 1
							break

				#octant 4
				elif initx>0 and inity<0 and initz>0:
					if len(poslist4) == 0:
						check3 = False; check1 = False; check8 = False;
						if initx < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check3) or np.any(check1) or np.any(check8)):
							poslist.append([initx, inity, initz]);
							poslist4.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist4_array = np.asarray(poslist4);
						dist = distance.cdist(poslist4_array, point);
						check = np.less(dist, 1); check3 = False; check1 = False; check8 = False;
						if initx < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if initz < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check3) or np.any(check1) or np.any(check8)):
							poslist.append([initx, inity, initz]);
							poslist4.append([initx, inity, initz]);
							index += 1
							break

				#octant 5
				elif initx>0 and inity>0 and initz<0:
					if len(poslist5) == 0:
						check6 = False; check8 = False; check1 = False;
						if initx < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check6) or np.any(check8) or np.any(check1)):
							poslist.append([initx, inity, initz]);
							poslist5.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist5_array = np.asarray(poslist5);
						dist = distance.cdist(poslist5_array, point);
						check = np.less(dist, 1); check6 = False; check8 = False; check1 = False;
						if initx < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check6) or np.any(check8) or np.any(check1)):
							poslist.append([initx, inity, initz]);
							poslist5.append([initx, inity, initz]);
							index += 1
							break

				#octant 6
				elif initx<0 and inity>0 and initz<0:
					if len(poslist6) == 0:
						check5 = False; check7 = False; check2 = False;
						if initx > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check5) or np.any(check7) or np.any(check2)):
							poslist.append([initx, inity, initz]);
							poslist6.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist6_array = np.asarray(poslist6);
						dist = distance.cdist(poslist6_array, point);
						check = np.less(dist, 1); check5 = False; check7 = False; check2 = False;
						if initx > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if inity < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check5) or np.any(check7) or np.any(check2)):
							poslist.append([initx, inity, initz]);
							poslist6.append([initx, inity, initz]);
							index += 1
							break

				#octant 7
				elif initx<0 and inity<0 and initz<0:
					if len(poslist7) == 0:
						check8 = False; check6 = False; check3 = False;
						if initx > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check8) or np.any(check6) or np.any(check3)):
							poslist.append([initx, inity, initz]);
							poslist7.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist7_array = np.asarray(poslist7);
						dist = distance.cdist(poslist7_array, point);
						check = np.less(dist, 1); check8 = False; check6 = False; check3 = False;
						if initx > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check8) or np.any(check6) or np.any(check3)):
							poslist.append([initx, inity, initz]);
							poslist7.append([initx, inity, initz]);
							index += 1
							break

				#octant 8
				elif initx>0 and inity<0 and initz<0:
					if len(poslist8) == 0:
						check7 = False; check5 = False; check4 = False;
						if initx < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check7) or np.any(check5) or np.any(check4)):
							poslist.append([initx, inity, initz]);
							poslist8.append([initx, inity, initz]);
							index += 1
							break
					else:
						poslist8_array = np.asarray(poslist8);
						dist = distance.cdist(poslist8_array, point);
						check = np.less(dist, 1); check7 = False; check5 = False; check4 = False;
						if initx < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if inity > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if initz > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check7) or np.any(check5) or np.any(check4)):
							poslist.append([initx, inity, initz]);
							poslist8.append([initx, inity, initz]);
							index += 1
							break


		#build up the rest of the protein
		for i in range(pro_len-1):

			while True:

				#candidate position of new particle
				pi = np.random.rand(1) * 2*np.pi;
				z_it = np.random.uniform(low=-1, high=1);
				x = np.sqrt(1-(z_it**2))*np.cos(pi) + poslist[index][0];
				y = np.sqrt(1-(z_it**2))*np.sin(pi) + poslist[index][1];
				z = z_it + poslist[index][2];
				point = np.array([[x, y, z]]);

				#check which octant candidate particle occupies and check for overlaps in corresponding list
				#octant 1
				if x>0 and y>0 and z>0:
					if len(poslist1) == 0:
						check2 = False; check4 = False; check5 = False;
						if x < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check2) or np.any(check4) or np.any(check5)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist1.append([x, y, z]);
								index += 1
								break
					else:
						poslist1_array = np.asarray(poslist1);
						dist = distance.cdist(poslist1_array, point);
						check = np.less(dist, 1); check2 = False; check4 = False; check5 = False;
						if x < 0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check2) or np.any(check4) or np.any(check5)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist1.append([x, y, z]);
								index += 1
								break

				#octant 2
				elif x<0 and y>0 and z>0:
					if len(poslist2) == 0:
						check1 = False; check3 = False; check6 = False;
						if x > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check1) or np.any(check3) or np.any(check6)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist2.append([x, y, z]);
								index += 1
								break
					else:
						poslist2_array = np.asarray(poslist2);
						dist = distance.cdist(poslist2_array, point);
						check = np.less(dist, 1); check1 = False; check3 = False; check6 = False;
						if x > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check1) or np.any(check3) or np.any(check6)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist2.append([x, y, z]);
								index += 1
								break

				#octant 3
				elif x<0 and y<0 and z>0:
					if len(poslist3) == 0:
						check4 = False; check2 = False; check7 = False;
						if x > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check4) or np.any(check2) or np.any(check7)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist3.append([x, y, z]);
								index += 1
								break
					else:
						poslist3_array = np.asarray(poslist3);
						dist = distance.cdist(poslist3_array, point);
						check = np.less(dist, 1); check4 = False; check2 = False; check7 = False;
						if x > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check4) or np.any(check2) or np.any(check7)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist3.append([x, y, z]);
								index += 1
								break

				#octant 4
				elif x>0 and y<0 and z>0:
					if len(poslist4) == 0:
						check3 = False; check1 = False; check8 = False;
						if x < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check3) or np.any(check1) or np.any(check8)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist4.append([x, y, z]);
								index += 1
								break
					else:
						poslist4_array = np.asarray(poslist4);
						dist = distance.cdist(poslist4_array, point);
						check = np.less(dist, 1); check3 = False; check1 = False; check8 = False;
						if x < 0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if z < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check3) or np.any(check1) or np.any(check8)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist4.append([x, y, z]);
								index += 1
								break

				#octant 5
				elif x>0 and y>0 and z<0:
					if len(poslist5) == 0:
						check6 = False; check8 = False; check1 = False;
						if x < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check6) or np.any(check8) or np.any(check1)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist5.append([x, y, z]);
								index += 1
								break
					else:
						poslist5_array = np.asarray(poslist5);
						dist = distance.cdist(poslist5_array, point);
						check = np.less(dist, 1); check6 = False; check8 = False; check1 = False;
						if x < 0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist1) > 0:
								poslist1_array = np.asarray(poslist1);
								dist = distance.cdist(poslist1_array, point);
								check1 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check6) or np.any(check8) or np.any(check1)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist5.append([x, y, z]);
								index += 1
								break

				#octant 6
				elif x<0 and y>0 and z<0:
					if len(poslist6) == 0:
						check5 = False; check7 = False; check2 = False;
						if x > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check5) or np.any(check7) or np.any(check2)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist6.append([x, y, z]);
								index += 1
								break
					else:
						poslist6_array = np.asarray(poslist6);
						dist = distance.cdist(poslist6_array, point);
						check = np.less(dist, 1); check5 = False; check7 = False; check2 = False;
						if x > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if y < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist2) > 0:
								poslist2_array = np.asarray(poslist2);
								dist = distance.cdist(poslist2_array, point);
								check2 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check5) or np.any(check7) or np.any(check2)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist6.append([x, y, z]);
								index += 1
								break

				#octant 7
				elif x<0 and y<0 and z<0:
					if len(poslist7) == 0:
						check8 = False; check6 = False; check3 = False;
						if x > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check8) or np.any(check6) or np.any(check3)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist7.append([x, y, z]);
								index += 1
								break
					else:
						poslist7_array = np.asarray(poslist7);
						dist = distance.cdist(poslist7_array, point);
						check = np.less(dist, 1); check8 = False; check6 = False; check3 = False;
						if x > -0.5:
							if len(poslist8) > 0:
								poslist8_array = np.asarray(poslist8);
								dist = distance.cdist(poslist8_array, point);
								check8 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist6) > 0:
								poslist6_array = np.asarray(poslist6);
								dist = distance.cdist(poslist6_array, point);
								check6 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist3) > 0:
								poslist3_array = np.asarray(poslist3);
								dist = distance.cdist(poslist3_array, point);
								check3 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check8) or np.any(check6) or np.any(check3)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist7.append([x, y, z]);
								index += 1
								break

				#octant 8
				elif x>0 and y<0 and z<0:
					if len(poslist8) == 0:
						check7 = False; check5 = False; check4 = False;
						if x < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check7) or np.any(check5) or np.any(check4)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist8.append([x, y, z]);
								index += 1
								break
					else:
						poslist8_array = np.asarray(poslist8);
						dist = distance.cdist(poslist8_array, point);
						check = np.less(dist, 1); check7 = False; check5 = False; check4 = False;
						if x < 0.5:
							if len(poslist7) > 0:
								poslist7_array = np.asarray(poslist7);
								dist = distance.cdist(poslist7_array, point);
								check7 = np.less(dist, 1);
						if y > -0.5:
							if len(poslist5) > 0:
								poslist5_array = np.asarray(poslist5);
								dist = distance.cdist(poslist5_array, point);
								check5 = np.less(dist, 1);
						if z > -0.5:
							if len(poslist4) > 0:
								poslist4_array = np.asarray(poslist4);
								dist = distance.cdist(poslist4_array, point);
								check4 = np.less(dist, 1);
						if ~(np.any(check) or np.any(check7) or np.any(check5) or np.any(check4)):
							if np.absolute(x)<extent/2-.5 and np.absolute(y)<extent/2-.5 and np.absolute(z)<extent/2-.5:
								poslist.append([x, y, z]);
								poslist8.append([x, y, z]);
								index += 1
								break




	return poslist



