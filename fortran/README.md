Trimmed version of SVN:merged_stilt_hysplit_v3 from Derek Mallia

===============================================================================================

THIS IS VERSION 2 OF THE STILT MODEL. INCLUDES FIXES FOR THE FOLLOWING SUBROUTINES:

adv2nt.f                     
adv3nt.f     
adv3ntZML.f   



WHAT IS FIXED: (FIX STARTING IN VERSION 1)

   Need to make corrections to several of the fortran subroutines in order to
   have the STILT model run correctly. Add the following line to the listed
   FORTRAN programs: adv2nt.f  adv3nt.f adv3ntZML.f
   
   	I1P=min(I1+1,size(s,1))
   	J1P=min(J1+1,size(s,2))
	
	
   These lines need to be added to each of the listed subroutinesright after 
   the commented line "!set upper index" and before the commented line 
   "!cyclinic boundary conditions". Then comment the following lines out with 
   a leading "!":
   
   	!I1P=I1+1
   	!J1P=J1+1
	


NEW FIXES:


     Additional set for STILT users who plan on using WRF-ARW data for the metfields. There is a bug in the STILT model
    that will cause particles to randomly terminate with high-resolution WRF met fields. To prevent this issue 
    Thomas Nehrkorn has modified the advrng.f subroutine within the STILT source code (this is located in your
    trunk/merged_stilt_hysplit directory. To edit this file go to the following directory copy the NEW advrng.f 
    with your old one. 

    	/uufs/chpc.utah.edu/common/home/lin-group1/model_inventory/WRF-STILT_patch

	Copy....

        cp advrng_new.f <your path>/trunk/merged_stilt_hysplit

    Recompile STILT when done. 

    After this retrieve the new stiltR trajec.r code (from the WRF-STILT_patch directory) which adds a new variable called mgmin. 
    Set this to >1000 and set your delta variable =< 5 . Copy this to your stiltR and edit your run script accordingly
    directory (i.e add mgmin variable to the Trajec.r function that you call. 

    This should mostly eliminate particle termination errors. For reference for WRF-STILT sample script in the 
    WRF-STILT_patch directory.







UPDATED BY DVM ON 10/12/13


