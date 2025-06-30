#!/bin/bash/
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

szr16jdjsgpd:q5=v0jJ

 out=Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+'GE_e4_merge_AGNseed001_SimBKG_CLUseed001'+'.fits'))
 ero_de = (out['OWNER']==2)|(out['OWNER']==0)
 to_process = ((out['OWNER'] == 2) | (out['OWNER'] == 0)) & (out['has_merged_events']) & (out['has_Sc1Cat'] == False)
 already_done = ((out['OWNER'] == 2) | (out['OWNER'] == 0)) & (out['has_merged_events']) & (out['has_Sc1Cat'])
 todo = (ero_de)&(to_process==False)&(already_done==False)
 out['fail'] = 0
 out['fail'][todo] = np.array(fails)
 out.write(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+'GE_e4_merge_AGNseed001_SimBKG_CLUseed001'+'_withFailReason.fits'), overwrite = True)
 out = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+'GE_e4_merge_AGNseed001_SimBKG_CLUseed001'+'_withFailReason.fits'))
 tosync_srv_map = out['SRVMAP'][(out['fail']==1)]
 for el in tosync_srv_map:
  	str_field = str(el).zfill(6)
  	print('mv ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'c030', '*') + ' ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_c030') )
 for el in tosync_srv_map:
  	str_field = str(el).zfill(6)
  	print('mv ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'eSASS', '*') + ' ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_eSASS') )

cd ~/workspace/erosim/Uchuu/LCerass
rsync -avz joco@raven.mpcdf.mpg.de:~/ptmp_joco/mpecl/comparat/data_s4_c030/00???? . # DONE
rsync -avz joco@raven.mpcdf.mpg.de:~/ptmp_joco/mpecl/comparat/data_s4_c030/0????? . # TODO

mv /home/idies/workspace/erosim/Uchuu/LCerass/096054/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096054/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090057/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090057/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094057/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094057/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097057/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097057/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084060/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/084060/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087060/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/087060/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091060/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/091060/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094060/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094060/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098060/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098060/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078063/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078063/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082063/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/082063/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085063/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085063/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088063/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/088063/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092063/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092063/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095063/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095063/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098063/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098063/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/076066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/079066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099066/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/099066/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/084069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/088069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/091069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097069/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097069/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/070072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099072/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/099072/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/069075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/082075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/088075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/091075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097075/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097075/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/069078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/084078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/087078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099078/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/099078/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098081/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098081/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098084/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098084/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098087/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098087/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098090/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098090/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098093/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098093/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098096/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098096/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/038099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098099/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098099/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/032102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/038102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/069102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/084102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/087102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099102/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/099102/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/029105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/029105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/032105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/038105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/042105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/048105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/051105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/054105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/057105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/069105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/082105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/088105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/091105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097105/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097105/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/023108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/027108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/030108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/033108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/036108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/036108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/039108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/042108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/049108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/052108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/055108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/058108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/061108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/061108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/064108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/064108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/070108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099108/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/099108/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/021111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/024111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/027111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/030111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/033111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/037111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/040111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/040111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/043111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/046111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/046111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/049111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/084111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/088111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/091111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097111/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097111/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/015114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/018114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/018114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/021114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/024114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/028114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/031114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/034114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/037114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/054114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/057114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/070114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/073114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/076114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/079114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099114/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/099114/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/008117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/008117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/012117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/015117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/018117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/018117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/022117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/025117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/028117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/032117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/038117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/042117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/048117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/052117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/055117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/058117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/082117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/088117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098117/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098117/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/005120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/009120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/012120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/015120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/019120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/022120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/026120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/029120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/029120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/033120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/036120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/036120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/039120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/043120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/046120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/046120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/057120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/070120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/084120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/087120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/091120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098120/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098120/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/005123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/009123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/012123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/016123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/019123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/023123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/026123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/030123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/034123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/037123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/048123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/051123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/055123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/058123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/069123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/076123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/079123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097123/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097123/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/005126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/009126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/013126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/013126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/016126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/020126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/024126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/027126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/031126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/038126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/042126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/049126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/064126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/064126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/075126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/082126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096126/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096126/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/006129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/009129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/013129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/013129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/017129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/017129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/021129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/025129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/028129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/032129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/036129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/036129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/040129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/040129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/051129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/055129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/070129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/078129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097129/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097129/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/006132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/010132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/010132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/014132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/018132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/018132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/022132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/026132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/030132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/034132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/038132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/042132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/049132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/057132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/061132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/061132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/065132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/069132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/073132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097132/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097132/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/006135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/010135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/010135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/014135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/019135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/023135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/027135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/031135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/039135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/043135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/048135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/052135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/064135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/064135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097135/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/097135/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/007138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/007138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/011138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/011138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/015138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/020138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/024138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/028138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/033138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/037138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/046138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/046138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/054138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/063138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/076138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098138/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098138/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/007141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/007141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/012141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/016141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/021141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/025141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/030141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/039141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/044141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/048141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/053141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/058141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/076141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099141/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/099141/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/002144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/007144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/007144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/012144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/017144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/017144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/022144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/027144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/032144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/037144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/042144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/052144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/057144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/062144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/072144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/076144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/081144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/091144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096144/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096144/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/003147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/008147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/008147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/013147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/013147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/019147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/024147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/029147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/029147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/034147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/040147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/040147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/061147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/061147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/082147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/087147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098147/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098147/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/003150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/009150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/014150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/020150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/026150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/031150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/037150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/043150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/049150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/054150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/083150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/089150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094150/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/094150/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/003153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/009153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/016153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/022153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/028153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/047153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/054153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/073153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/079153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/085153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/092153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098153/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098153/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/003156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/010156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/010156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/017156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/017156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/024156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/031156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/038156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/052156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/066156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/073156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/080156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/087156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093156/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/093156/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/004159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/004159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/012159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/020159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/027159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/035159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/043159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/051159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/067159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/082159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098159/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/098159/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/005162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/014162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/023162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/032162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/050162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/059162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/068162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/077162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/086162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095162/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/095162/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/005165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/016165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/026165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/037165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/048165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/058165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/069165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/079165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090165/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090165/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/006168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/019168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/032168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/045168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/058168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/071168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/084168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096168/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/096168/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/008171/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/008171/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025171/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/025171/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041171/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/041171/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057171/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/057171/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074171/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/074171/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090171/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/090171/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/011174/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/011174/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034174/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/034174/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056174/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/056174/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079174/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/079174/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020177/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/020177/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060177/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/060177/s4_c030/
mv /home/idies/workspace/erosim/Uchuu/LCerass/001180/c030/* /home/idies/workspace/erosim/Uchuu/LCerass/001180/s4_c030/


mv /home/idies/workspace/erosim/Uchuu/LCerass/096054/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096054/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090057/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090057/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094057/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094057/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097057/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097057/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084060/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/084060/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087060/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/087060/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091060/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/091060/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094060/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094060/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098060/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098060/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078063/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078063/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082063/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/082063/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085063/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085063/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088063/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/088063/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092063/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092063/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095063/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095063/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098063/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098063/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/076066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/079066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099066/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/099066/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/084069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/088069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/091069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097069/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097069/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/070072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099072/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/099072/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/069075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/082075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/088075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/091075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097075/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097075/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/069078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/084078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/087078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099078/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/099078/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098081/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098081/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098084/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098084/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098087/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098087/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098090/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098090/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098093/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098093/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098096/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098096/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/038099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098099/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098099/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/032102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/038102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/069102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/084102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/087102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099102/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/099102/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/029105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/029105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/032105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/038105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/042105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/048105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/051105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/054105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/057105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/069105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/082105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/088105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/091105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097105/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097105/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/023108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/027108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/030108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/033108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/036108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/036108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/039108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/042108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/049108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/052108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/055108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/058108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/061108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/061108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/064108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/064108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/070108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099108/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/099108/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/021111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/024111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/027111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/030111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/033111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/037111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/040111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/040111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/043111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/046111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/046111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/049111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/084111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/088111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/091111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097111/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097111/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/015114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/018114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/018114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/021114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/024114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/028114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/031114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/034114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/037114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/054114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/057114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/070114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/073114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/076114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/079114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099114/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/099114/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/008117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/008117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/012117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/015117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/018117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/018117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/022117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/025117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/028117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/032117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/038117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/042117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/048117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/052117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/055117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/058117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/082117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/088117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/088117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098117/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098117/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/005120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/009120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/012120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/015120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/019120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/022120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/026120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/029120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/029120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/033120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/036120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/036120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/039120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/043120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/046120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/046120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/057120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/070120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/084120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/087120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/091120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098120/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098120/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/005123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/009123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/012123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/016123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/019123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/023123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/026123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/030123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/034123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/037123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/048123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/051123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/055123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/058123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/069123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/076123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/079123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097123/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097123/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/005126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/009126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/013126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/013126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/016126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/020126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/024126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/027126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/031126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/038126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/042126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/049126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/064126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/064126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/075126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/075126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/082126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096126/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096126/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/006129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/009129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/013129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/013129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/017129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/017129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/021129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/025129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/028129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/032129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/036129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/036129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/040129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/040129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/051129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/055129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/055129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/070129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/070129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/078129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/078129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097129/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097129/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/006132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/010132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/010132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/014132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/018132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/018132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/022132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/026132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/030132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/034132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/038132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/042132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/049132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/057132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/061132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/061132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/065132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/065132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/069132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/073132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097132/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097132/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/006135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/010135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/010135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/014135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/019135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/023135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/027135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/031135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/039135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/043135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/048135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/052135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/064135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/064135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/097135/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/097135/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/007138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/007138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/011138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/011138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/015138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/015138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/020138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/024138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/028138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/033138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/033138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/037138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/046138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/046138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/054138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/063138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/063138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/076138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098138/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098138/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/007141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/007141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/012141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/016141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/021141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/021141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/025141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/030141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/030141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/039141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/039141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/044141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/044141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/048141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/053141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/053141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/058141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/076141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/099141/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/099141/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/002144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/002144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/007144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/007144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/012144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/017144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/017144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/022144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/027144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/032144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/037144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/042144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/042144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/052144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/057144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/062144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/062144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/072144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/072144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/076144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/076144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/081144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/081144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/091144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/091144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096144/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096144/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/003147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/008147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/008147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/013147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/013147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/019147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/024147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/029147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/029147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/034147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/040147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/040147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/061147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/061147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/082147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/087147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098147/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098147/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/003150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/009150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/014150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/020150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/026150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/031150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/037150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/043150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/049150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/049150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/054150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/083150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/083150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/089150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/089150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/094150/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/094150/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/003153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/009153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/009153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/016153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/022153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/022153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/028153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/028153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/047153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/047153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/054153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/054153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/073153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/079153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/085153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/085153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/092153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/092153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098153/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098153/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/003156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/003156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/010156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/010156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/017156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/017156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/024156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/024156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/031156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/031156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/038156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/038156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/052156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/052156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/066156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/066156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/073156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/073156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/080156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/080156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/087156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/087156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/093156/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/093156/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/004159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/004159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/012159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/012159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/020159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/027159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/027159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/035159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/035159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/043159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/043159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/051159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/051159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/067159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/067159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/082159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/082159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/098159/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/098159/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/005162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/014162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/014162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/023162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/023162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/032162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/050162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/050162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/059162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/059162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/068162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/068162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/077162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/077162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/086162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/086162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/095162/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/095162/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/005165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/005165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/016165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/016165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/026165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/026165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/037165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/037165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/048165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/048165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/058165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/069165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/069165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/079165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090165/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090165/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/006168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/006168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/019168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/019168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/032168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/032168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/045168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/045168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/058168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/058168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/071168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/071168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/084168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/084168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/096168/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/096168/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/008171/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/008171/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/025171/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/025171/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/041171/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/041171/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/057171/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/057171/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/074171/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/074171/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/090171/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/090171/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/011174/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/011174/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/034174/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/034174/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/056174/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/056174/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/079174/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/079174/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/020177/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/020177/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/060177/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/060177/s4_eSASS/
mv /home/idies/workspace/erosim/Uchuu/LCerass/001180/eSASS/* /home/idies/workspace/erosim/Uchuu/LCerass/001180/s4_eSASS/

BG model files are here :
~/workspace/erosim/simput/bkg_erosita_simput_full_sky
catalog.fits has the full map

# make a sky map with input files
python make_summarySimEvt_skymap.py # DONE, OK, all files are there.

# in any container with python
nohup python merge_events_onlyBG.py  > logs/merge_events_onlyBG.log  & # ONGOING

nohup python merge_events.py 1 1     > logs/merge_events_1_1.log  & # DONE
nohup python merge_events_noCLU.py 1 > logs/merge_events_noCLU_1.log & # TODO
nohup python merge_events_noAGN.py 1 > logs/merge_events_noAGN_1.log & # TODO

nohup python merge_events.py 2 2     > logs/merge_events_2_2.log  & # DONE
nohup python merge_events_noCLU.py 2 > logs/merge_events_noCLU_2.log & # TODO
nohup python merge_events_noAGN.py 2 > logs/merge_events_noAGN_2.log & # TODO

nohup python merge_events.py 3 3     > logs/merge_events_3_3.log  & # ONGOING
nohup python merge_events_noCLU.py 3 > logs/merge_events_noCLU_3.log & # TODO
nohup python merge_events_noAGN.py 3 > logs/merge_events_noAGN_3.log & # TODO

nohup python merge_events.py 4 4     > logs/merge_events_4_4.log  & # ONGOING
nohup python merge_events_noCLU.py 4 > logs/merge_events_noCLU_4.log & # TODO
nohup python merge_events_noAGN.py 4 > logs/merge_events_noAGN_4.log & # TODO

nohup python merge_events.py 5 5     > logs/merge_events_5_5.log  & # ONGOING
nohup python merge_events_noCLU.py 5 > logs/merge_events_noCLU_5.log & # TODO
nohup python merge_events_noAGN.py 5 > logs/merge_events_noAGN_5.log & # TODO

nohup python merge_events.py 6 6     > logs/merge_events_6_6.log  & # ONGOING
nohup python merge_events_noCLU.py 6 > logs/merge_events_noCLU_6.log & # TODO
nohup python merge_events_noAGN.py 6 > logs/merge_events_noAGN_6.log & # TODO

nohup python merge_events.py 7 7     > logs/merge_events_7_7.log  & # ONGOING
nohup python merge_events_noCLU.py 7 > logs/merge_events_noCLU_7.log & # TODO
nohup python merge_events_noAGN.py 7 > logs/merge_events_noAGN_7.log & # TODO

nohup python merge_events.py 8 8     > logs/merge_events_8_8.log  & # ONGOING
nohup python merge_events_noCLU.py 8 > logs/merge_events_noCLU_8.log & # TODO
nohup python merge_events_noAGN.py 8 > logs/merge_events_noAGN_8.log & # TODO


GE_e4_merge_AGNseed00?_SimBKG : AGN + BKG
GE_e4_merge_AGNseed00?_SimBKG_CLUseed00? : AGN + BKG + Cluster
GE_e4_merge_SimBKG : only BKG
GE_e4_merge_SimBKG_CLUseed00? : BKG + Cluster

# write eSASS commands in the relevant folders
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG                       > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG.log                       & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed001            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed001.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG_CLUseed002.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed002            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed002.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG_CLUseed003.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed003            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed003.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed004_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed004_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed004_SimBKG_CLUseed004.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed004            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed004.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed005_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed005_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed005_SimBKG_CLUseed005.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed005            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed005.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed006_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed006_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed006_SimBKG_CLUseed006.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed006            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed006.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed007_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed007_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed007_SimBKG_CLUseed007.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed007            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed007.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed008_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed008_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed008_SimBKG_CLUseed008.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed008            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed008.log            & # DONE


#monitor
# monitor progress
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

nohup python make_summary_skymap.py GE_e4_merge_SimBKG                       > logs/summary_sky_map_GE_e4_merge_SimBKG.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed001_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed001_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/summary_sky_map_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed001            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed001.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed002_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed002_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 > logs/summary_sky_map_GE_e4_merge_AGNseed002_SimBKG_CLUseed002.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed002            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed002.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed003_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed003_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 > logs/summary_sky_map_GE_e4_merge_AGNseed003_SimBKG_CLUseed003.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed003            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed003.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed004_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed004_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 > logs/summary_skymap_GE_e4_merge_AGNseed004_SimBKG_CLUseed004.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed004            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed004.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed005_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed005_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 > logs/summary_skymap_GE_e4_merge_AGNseed005_SimBKG_CLUseed005.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed005            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed005.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed006_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed006_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 > logs/summary_skymap_GE_e4_merge_AGNseed006_SimBKG_CLUseed006.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed006            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed006.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed007_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed007_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 > logs/summary_skymap_GE_e4_merge_AGNseed007_SimBKG_CLUseed007.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed007            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed007.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed008_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed008_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 > logs/summary_skymap_GE_e4_merge_AGNseed008_SimBKG_CLUseed008.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed008            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed008.log            & # TODO


# write all commands for the list of folders above.
python write_exec_loop.py # > exec.sh # TODO

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/runs

# execute all commands of interest in an eSASS loaded container

nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0000.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0050.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0100.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0150.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0200.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0250.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0300.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0350.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0400.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0450.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0500.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0550.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0600.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0650.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0650.log & # ONGOING

nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0000.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0050.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0100.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0150.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0200.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0250.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0300.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0350.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0400.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0450.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0500.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0550.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0600.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0650.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0700.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0700.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0750.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0750.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0800.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0800.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0850.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0850.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0900.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0900.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0950.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0950.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1000.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1050.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1100.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1150.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1200.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1250.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1300.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1350.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1400.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1450.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1450.log & # ONGOING

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/runs

nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0000.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0050.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0100.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0150.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0200.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0250.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0300.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0350.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0400.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0450.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0500.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0550.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0600.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0650.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0700.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0700.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0750.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0750.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0800.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0800.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0850.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0850.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0900.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0900.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0950.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0950.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1000.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1050.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1100.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1150.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1200.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1250.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1300.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1350.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1400.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1450.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1450.log & # ONGOING


export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/runs

nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0000.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0050.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0100.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0150.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0200.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0250.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0300.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0350.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0400.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0450.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0500.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0550.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0600.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0650.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0700.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0700.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0750.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0750.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0800.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0800.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0850.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0850.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0900.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0900.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0950.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0950.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1000.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1050.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1100.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1150.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1200.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1250.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1300.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1350.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1400.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1450.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1450.log & # ONGOING


nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0000.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0050.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0100.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0150.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0200.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0250.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0300.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0350.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0400.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0450.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0500.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0550.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0600.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0650.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0700.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0700.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0750.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0750.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0800.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0800.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0850.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0850.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0900.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0900.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_0950.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_0950.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1000.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1050.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1100.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1150.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1200.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1250.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1300.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1350.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1400.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1450.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1500.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1550.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1600.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1650.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_processing_1700.sh > logs/GE_e4_merge_AGNseed002_SimBKG_processing_1700.log & # ONGOING

export UCHUU='/home/idies/workspace/erosim/Uchuu'                         
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'           
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data' 
cd $GIT_STMOD/src/esass/runs                                              


nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0000.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0000.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0050.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0050.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0100.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0100.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0150.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0150.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0200.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0200.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0250.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0250.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0300.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0300.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0350.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0350.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0400.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0400.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0450.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0450.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0500.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0500.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0550.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0550.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0600.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0600.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0650.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0650.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0700.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0700.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0750.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0750.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0800.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0800.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0850.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0850.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0900.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0900.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_0950.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_0950.log &# ONGOING

nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1000.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1000.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1050.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1050.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1100.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1100.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1150.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1150.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1200.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1200.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1250.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1250.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1300.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1300.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1350.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1350.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1400.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1400.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1450.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1450.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1500.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1500.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1550.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1550.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1600.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1600.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1650.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1650.log &# ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_processing_1700.sh > logs/GE_e4_merge_AGNseed003_SimBKG_processing_1700.log &# ONGOING

# TODO
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/runs

nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0000.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0000.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0050.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0050.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0100.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0100.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0150.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0150.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0200.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0200.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0250.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0250.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0300.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0300.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0350.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0350.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0400.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0400.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0450.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0450.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0500.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0500.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0550.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0550.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0600.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0600.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0650.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0650.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0700.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0700.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0750.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0750.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0800.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0800.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0850.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0850.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0900.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0900.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_0950.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_0950.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1000.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1000.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1050.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1050.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1100.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1100.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1150.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1150.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1200.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1200.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1250.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1250.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1300.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1300.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1350.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1350.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1400.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1400.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1450.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1450.log &
nohup sh GE_e4_merge_SimBKG_CLUseed003_processing_1500.sh > logs/GE_e4_merge_SimBKG_CLUseed003_processing_1500.log &

# TODO
# start realizations 5,6,7,8
# other scaling relation sim


export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/


python create_summary_files_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed001
python create_summary_files_RS.py GE_e4_merge_AGNseed001_SimBKG
python create_summary_files_RS.py GE_e4_merge_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed002
python create_summary_files_RS.py GE_e4_merge_AGNseed002_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed003
python create_summary_files_RS.py GE_e4_merge_AGNseed003_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed004
python create_summary_files_RS.py GE_e4_merge_AGNseed004_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed005
python create_summary_files_RS.py GE_e4_merge_AGNseed005_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed006
python create_summary_files_RS.py GE_e4_merge_AGNseed006_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed007
python create_summary_files_RS.py GE_e4_merge_AGNseed007_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed008
python create_summary_files_RS.py GE_e4_merge_AGNseed008_SimBKG

python summary_plotting_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins

python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed001 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed002 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed003 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed004 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed005 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed006 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed007 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed008 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed001_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed002_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed003_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed004_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed005_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed006_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed007_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed008_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins














#
# test fields and scripts :
#
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
cd $GIT_STMOD/src/esass/runs
sh GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.sh  > logs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.log & # running
sh GE_e4_merge_SimBKG_CLUseed001_processing_0000.sh             > logs/GE_e4_merge_SimBKG_CLUseed001_processing_0000.log            & # running
sh GE_e4_merge_SimBKG_processing_0000.sh                        > logs/GE_e4_merge_SimBKG_processing_0000.log                       & # running


cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd /home/idies/workspace/erosim/software/st_mod/src/esass
python photon_matching_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 121048


# test field
# execute eSASS commands on test field:
cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd $GIT_STMOD/src/esass
python photon_matching_RS.py GE_e4_merge_AGNseed001_SimBKG 121048

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd $GIT_STMOD/src/esass
python photon_matching_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 121048

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_SimBKG/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
#python photon_matching_RS.py GE_e4_merge_SimBKG 121048

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd $GIT_STMOD/src/esass
python photon_matching_RS.py GE_e4_merge_SimBKG_CLUseed001 121048

