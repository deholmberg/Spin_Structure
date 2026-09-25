
#ifndef BACKGROUND_RUNS_H
#define BACKGROUND_RUNS_H

    vector<int> Su22_F_Runs = {16194}; // Foils only
    vector<int> Su22_C_Runs = {16290, 16291, 16292, 16293, 16296, 16297}; // C runs
    vector<int> Su22_CH2_Runs={16298, 16299, 16300, 16302, 16303}; // CH2 runs
    vector<int> Su22_ET_Runs ={16309}; // ET runs
    //16700, 16701, 16702, 16704 // More C runs
    
    // Stable runs for the negative solenoid configuration of Fa22
    vector<int> Fa22Neg_ET_Runs = {16975}; // Empty target run 
    vector<int> Fa22Neg_F_Runs =  {16979}; // Foils only run
    vector<int> Fa22Neg_C_Runs =  {17180, 17181, 17182, 17183}; // C runs
    vector<int> Fa22Neg_CH2_Runs ={17132, 17133, 17134, 17135}; // CH2 runs

    // Stable runs for the positive solenoid configuration of Fa22
    vector<int> Fa22Pos_C_Runs = {17316, 17317, 17319, 17320}; // C runs
    vector<int> Fa22Pos_CH2_Runs={17387, 17388, 17393}; // CH2 runs

    // Stable runs for the inbending torus configuration of Sp23Inb
    vector<int> Sp23Inb_F_Runs = {17763, 17764, 17765}; // Foil runs
    vector<int> Sp23Inb_ET_Runs= {17766, 17767, 17768}; // Empty target runs
    vector<int> Sp23Inb_C_Runs = {17566, 17567, 17569, 17744, 17745, 17746}; // C runs
    vector<int> Sp23Inb_CH2_Runs={17625, 17626, 17627, 17635, 17637, 17638}; // CH2 runs
	//17514, 17516, 17517, 17518, 17520, 
	//17668, 17669, 17670, 17671, 17672, 17673, 17674, 17675, 17676 
	//17516, 17517, 17520, 17669, 17670, 17672  // CD2 runs
    //vector<int> Sp23Inb_CD2_Runs={17516, 17517, 17520}; // CD2 runs
    vector<int> Sp23Inb_CD2_Runs={17669, 17670, 17672, 17516, 17517, 17520}; // CD2 runs

#endif
