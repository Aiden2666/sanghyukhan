function [A, resolution, slice_z] = parse_pQCT (full_name)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
%  8/30/06 Derek Lindsey
%
%  Read the pertinent header information from the pQCT CT scans
%
%  VoxelSize - Byte 12          double
%  Patient Name - Byte 1099,    41 character string, 1st byte is the # of
%                               characters used in the string
%  X-Matrix Size of Image - Byte 1529   2 byte unsigned
%  Y-Matrix Size of Image - Byte 1531   2 byte unsigned

%  8/30/06 Derek Lindsey
%
%  Read the header from the pQCT CT scans

%  Data types 

%    boolean     : 1 byte boolean (0=false; 1=true)
%    byte        : 1 byte unsigned
%    shortint    : 1 byte signed
%    word        : 2 byte unsigned
%    integer     : 2 byte signed
%    longint     : 4 byte signed

%    single      : 4 byte floating point
%    real        : 6 byte floating point
%    double      : 8 byte floating point

%    char        : 1 byte character
%    string      : 256 byte of characters first character defines used length
%    string[nn]  : nn+1 byte of characters first character defines used length

%    pointer     : 4 byte pointer
%    ^anytype    : is a pointer to a variable of type "anytype"

% Basic Structure of CT Image Files

% FilePrefix    (Header Prefix Type)
% DetRec        (Describes the geometry of detectors...)
% PatInfoRec    (Information about the patient)
% PicInfoRec    (Information about the ct-image)
% CtImageData   (CT image data)

 fid = fopen(full_name, 'r');  
% 
% % FilePrefix
     HeadVers = fread(fid, 1, 'int32');                  % Version of header:   i.e.:3
     HeadLen = fread(fid, 1, 'int32');                   % Length of header

% %DetRec
     RecSize = fread(fid, 1, 'int32')                   % = SizeOf (DetRecType)
     VoxelSize = fread(fid, 1, 'double')                 % Size of a voxel [mm]
     NumIntv = fread(fid, 1, 'uint16')                  % total number of intervals
     AirRadius = fread(fid, 1, 'double')                % radius to take airvalues
     NumSlices = fread(fid, 1, 'uint16')                 % total number of slices
     SliceStart = fread(fid, 1, 'double')                % Z-Pos of first slice [mm]
     SliceDist = fread(fid, 1, 'double')                 % DeltaZ between two slices [mm]
     SlicePause_NotUsed = fread(fid, 1, 'uint16');       % Seconds pause between Slices
%     NumBlocks = fread(fid, 1, 'uint16');                % Number of Blocks
%     BlockPause_NotUsed = fread(fid, 1, 'uint16');       % Seconds pause between Blocks
%     WaitSec = fread(fid, 1, 'uint16');                  % Seconds to wait after HV is switched on
%     %Input
%     DistSourceCenter = fread(fid, 1, 'double');         % Distance source to rotation center [mm]
%     OffSetCenterMm = fread(fid, 1, 'double');           % x-axis: Distance zero-pos to center [mm]
%     AddValueMm = fread(fid, 1, 'double');               % value added to start and end of measrange
%     SLVoxelSizeRef = fread(fid, 1, 'double') ;          % reference voxelsize to calculate SLFactor
%     SLFactorRef = fread(fid, 1, 'double');              % reference SLFactor to calculate SLFactor
%     %Calculated
%     SLFactor = fread(fid, 1, 'double');                 % SLFactor = SLFactorRef/SLVoxelSizeRef*VoxelSize
%     %Input
%     DeadTimeMySec = fread(fid, 1, 'double');            % {Deadtime of detectors in microseconds    }
%     BeamHard_SF1 = fread(fid, 1, 'double');             % {Parameter SF1 for beam hardening corr.   }
%     BeamHard_SF2 = fread(fid, 1, 'double');             % {Parameter SF2 for beam hardening corr.   }
%     BeamHard_SF3 = fread(fid, 1, 'double');             % {Parameter SF3 for beam hardening corr.   }
%     BeamHard_SF4 = fread(fid, 1, 'double');             % {Parameter SF4 for beam hardening corr.   }
%     Delta_Theta = fread(fid, 1, 'double');              % {Angle between two detectors              }
%     DetDir = fread(fid, 1, 'int16');                    % {+1: increase Angle with DetNo; -1: invers}
%     Theta_Det1= fread(fid, 1, 'double');                % {Angle Theta of physical Det1             }
%     FirstUsedDet = fread(fid, 1, 'uint16');             % {First Physikal Detector to be used       }
%     NumDet = fread(fid, 1, 'uint16');                   % {Number of Dets to be used,               }
%     ImpPerMmX= fread(fid, 1, 'double');                 % {Impulse of angular encoder per mm; x-axis}
%     ImpPerDeg= fread(fid, 1, 'double');                 % {Impulse of angular encoder per deg;y-axis}
%     ImpPerMmZ= fread(fid, 1, 'double');                 % {Impulse of angular encoder per mm; z-axis}
%     TimRefScale = fread(fid, 1, 'double');              % {Counts of Time reference per InterValTime}
%     InterValTime = fread(fid, 1, 'double');             % {Time per normalized interval [sec]       }
%     DetNumByte = fread(fid, 1, 'uint8');                % {Number of bytes per dataintv. from device}
%     TimRefNumByte = fread(fid, 1, 'uint8');             % {Number of bytes per timereference data   }
%     XRefNumByte = fread(fid, 1, 'uint8');               % {Number of bytes per postition reference  }	n
%     DetEprFromDetPhys = fread(fid, MaxNumDet, 'uchar'); %	: Array[1..MaxNumDet] of Byte;	% {input; order of det}
%     StartAngle 	= fread(fid, 1, 'double');              % {Offset of Angle                          }
%     OddBlockDir = fread(fid, 1, 'int16');               % {+1: odd blocks increase Angle with ScanNo}
%                                                         % {-1: odd blocks decrease Angle with ScanNo}
%     OddScanDir = fread(fid, 1, 'int16');                % {+1: odd scans increase Xpos with IntvNo  }
%                                                         % {-1: odd scans decrease Xpos with IntvNo  }
%     SystemDir = fread(fid, 1, 'int16');                 % {+1: right handed system                  }
%                                                         % {-1: left handed system                   }
%     MirrorDir = fread(fid, 1, 'int16');                 % {+1: do not invert picture                }
%                                                         % {-1: invert picture                       }
%     % {--- Calculated ---}
%     LastUsedDet = fread(fid, 1, 'uint16');              % {[1] Output; (= FirstUsedDet + NumDet - 1)}
%     NumBScans = fread(fid, 1, 'uint16');                % {[1] Output; number of translations per block}
%     NumBProj = fread(fid, 1, 'uint16');                 % {[1] Output; number of projections per block}
%     NumBProjUsed = fread(fid, 1, 'uint16');             % {[1] Output; number of used proj. per block}
%     NumBProjUnused = fread(fid, 1, 'uint16');           % {[1] Output; number ou unused proj. per block}
%     NumAllProjUsed = fread(fid, 1, 'uint16');           % {[1] Output; number of all used projections}
%     StepAngle = fread(fid, 1, 'double');                % {[deg] Output; NumDet * Delta_Theta}
%     Delta_Angle = fread(fid, 1, 'double');              % {[deg] Output; Delta_Theta / NumBlocks}
%     IntvNumByte = fread(fid, 1, 'uint8');               % {[1] Output;}
%     % {-------------------------}
%     ZPos_Start = fread(fid, 1, 'double');               % {[mm] Input;       Start of SV               }
%     ZPos_DiffRefMeas = fread(fid, 1, 'double');         % {[mm] Output&Input Dist. Ref. to Meas pos.   }
%     ZPos_MeasCT = fread(fid, 1, 'double');              % {[mm] Input;       Meas. position            }
%     ZPos_Step = fread(fid, 1, 'double');                % {[mm] Input;       Line Dist. of SV          }
%     ObjLen = fread(fid, 1, 'double');                   % {[mm] Input;       Length of object          }
%     Percent = fread(fid, 1, 'double');                  % {[%]  Input; Dist. Ref. to Meas pos [%ObjLen]}
% 
%     ZPos_NumSteps = fread(fid, 1, 'uint16');            % {[1]          Input;  Number of SV-Lines   }
%     MSpeed_MmPerSec = fread(fid, 1, 'double');          % {[mm/Sec]     Output; (not used)}
%     MTime_mSecPerIntv = fread(fid, 1, 'double');        % {[mSec/Intv]  Output; (not used)}
%     MTime_SecPerScan = fread(fid, 1, 'double');         % {[Sec/Scan]   Output; (not used)}
%     MTime_SecPerSlice = fread(fid, 1, 'double');        % {[Sec/Slice]  Output; (not used)}
%     RateMin_NeverUsed = fread(fid, 1, 'double');        % {[Counts/Sec] Input;  (not used)}		n
%     RateMax_NeverUsed = fread(fid, 1, 'double');        % {[Counts/Sec] Input;  (not used)} 		n
% 
%     SampleTimeMySec = fread(fid, 1, 'double');          % {[mySec] Input; Sample time of motor controller}
%     XUserDir_NeverUsed = fread(fid, 1, 'int16');        % {[1]          Input;  (not used)} 		n
%     XUserOffs_NeverUsed	= fread(fid, 1, 'double');      % {[mm]         Input;  (not used)} 		n
%     XPos_Min_NeverUsed = fread(fid, 1, 'double');       % {[mm]     Input; (not used)} 			n
%     XPos_Max_NeverUsed = fread(fid, 1, 'double');       % {[mm]     Input; (not used)} 			n
%     XV_Min_NeverUsed = fread(fid, 1, 'double');         % {[mm/Sec] Input; (not used)} 			n
%     XV_Max_NeverUsed = fread(fid, 1, 'double');         % {[mm/Sec] Input; (not used)} 			n
%     XAcc = fread(fid, 1, 'double');                     % {[mm/Sec²]Input; X-Axis: Accelleration}
%     XV_0 = fread(fid, 1, 'double');                     % {[mm/Sec] Input; Velo Meas. minimum (not used)}
%     XV_1 = fread(fid, 1, 'double');                     % {[mm/Sec] Input; V. Measurement}
%     XV_2 = fread(fid, 1, 'double');                     % {[mm/Sec] Input; V. Measurement move}
%     XV_S = fread(fid, 1, 'double');                     % {[mm/Sec] Input; V. Search}
%     XV_M = fread(fid, 1, 'double');                     % {[mm/Sec] Input; V. Normal Move
% 
%     YUserDir_NeverUsed	= fread(fid, 1, 'int16');       % {[1]      Input}				n
%     YUserOffs_NeverUsed	= fread(fid, 1, 'double');      % {[mm]     Input}				n
%     YPos_Min_NeverUsed = fread(fid, 1, 'double');       % {[mm]     Input}				n
%     YPos_Max_NeverUsed = fread(fid, 1, 'double');       % {[mm]     Input}				n
%     YV_Min_NeverUsed = fread(fid, 1, 'double');         % {[mm/Sec] Input}				n
%     YV_Max_NeverUsed = fread(fid, 1, 'double');         % {[mm/Sec] Input}				n
%     YAcc = fread(fid, 1, 'double');                     % {[mm/Sec²]Input}
%     YV_0 = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     YV_1 = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     YV_2 = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     YV_S = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     YV_M = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
% 
%     ZUserDir_NeverUsed = fread(fid, 1, 'int16');        % {[1]      Input}				n
%     ZUserOffs_NeverUsed = fread(fid, 1, 'double');      % {[mm]     Input}				n
%     ZPos_Min_NeverUsed = fread(fid, 1, 'double');       % {[mm]     Input}				n
%     ZPos_Max_NeverUsed = fread(fid, 1, 'double');       % {[mm]     Input}				n
%     ZV_Min_NeverUsed = fread(fid, 1, 'double');         % {[mm/Sec] Input}				n
%     ZV_Max_NeverUsed = fread(fid, 1, 'double');         % {[mm/Sec] Input}				n
%     ZAcc = fread(fid, 1, 'double');                     % {[mm/Sec²]Input	}
%     ZV_0 = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     ZV_1 = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     ZV_2 = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     ZV_S = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     ZV_M = fread(fid, 1, 'double');                     % {[mm/Sec] Input}
%     MeasInfo_a = fread(fid, 81, '*char');               % Array[1..4] of String80;
%     MeasInfo_b = fread(fid, 81, '*char');               % Array[1..4] of String80;
%     MeasInfo_c = fread(fid, 81, '*char');               % Array[1..4] of String80;
%     MeasInfo_d = fread(fid, 81, '*char');               % Array[1..4] of String80;
%     MeasDate = fread(fid, 1, 'int32')
%     MeasType = fread(fid, 1, 'int32');                  % {0: SV; 1: CT}
%     MeasState = fread(fid, 1, 'int32');                 % {not used}
%     MeasQuality_Free = fread(fid, 1, 'int8');           % {not used}					*
%     MeasQuality_Ln = fread(fid, 1, 'int8');             % {LnCheck}					*
%     MeasQuality_Offset = fread(fid, 1, 'int8');         % {OffsetCheck}					*
%     MeasQuality_Range = fread(fid, 1, 'int8');          % {RangeCheck}					*
%                                                         % {MeasQuality: positive = inside limits; 
%                                                         %               negative = outside limits}
%     IntvNoNumByte = fread(fid, 1, 'uint8');             % {Input}
%     ChkSumNumByte = fread(fid, 1, 'uint8');             % {Input}
%     HexFlg = fread(fid, 1, 'int8');                     % {Input}
%     DetPhysFromDetEpr = fread(fid, MaxNumDet, 'int8');  %: Array[1..MaxNumDet] of Byte;	% {OUTPUT}
%     ColScale = fread(fid, 1, 'float32');
%     ColOffs = fread(fid, 1, 'float32');
%     % {Reserve 	: Array[1..32] of Byte;} {until 28.11.96 P.W}
%     % {---------------------------}
%     % {since 28.11.96 P.W:} 
%     QaFlg = fread(fid, 1, 'int8');                      % { 1 Byte }
%     Sv_Comp_Delta = fread(fid, 1, 'float32');           % { 4 Byte }
%     DeviceTyp = fread(fid, 13, '*char');                % {13 Byte!}
%     VoxSizChrCT = fread(fid, 1, '*char');               % { 1 Byte }
% 
%     Info_TubeCurrentMa = fread(fid, 1, 'float32')       % { 4 Byte }{nominal TubeCurrent [mA]}		*
%     Info_TubePotentialKv = fread(fid, 1, 'float32')     % { 4 Byte }{nominal TubePotential [kV]}	*
%     GroupNo = fread(fid, 1, 'uint8');                   % { 1 Byte }{Number of Meas. Group}		*
%     SliceNoGrp = fread(fid, 1, 'uint8');                % { 1 Byte }{Slice Number of Group}		*
%     DistFlg	= fread(fid, 1, 'int8');                    % { 1 Byte }{Percent or Distance mode}		*
%     AsymFlg	= fread(fid, 1, 'int8');                    % { 1 Byte }{Slice symmetric or not}		*
%     Reserve	= fread(fid, 1, 'int8');                    % Array[1..1] of Byte;					%	*
% 
%     
% %PatInfoRec (information about the patient)
% 
%     RecSize_ = fread(fid, 1, 'int32');                  %{= SizeOf(PatInfoRecType) }
%     PatGender = fread(fid, 1, 'uint16');                %{Gender of patient}
%     PatEthnicGrp = fread(fid, 1, 'uint16');             %{Number of ethnic group of patient}
%     PatMeasNo = fread(fid, 1, 'uint16')                 %{Patients meas number }
%     PatNo = fread(fid, 1, 'int32')                      %{Patient number }
%     PatBirth = fread(fid, 1, 'int32')                   %{date of birth of patient}
%     PatMenoAge = fread(fid, 1, 'int32')                 %{meno age of patient}
%     PatName = fread(fid, 41, '*char')';                 %{Patient familyname & firstname }
%     PatNameChar = int8(PatName(1));                     %The first byte is the number of characters used
%     PatName = PatName(2:(1+PatNameChar))                %Pulling the Patient Name out
%     ZeroChar1 = fread(fid, 1, '*char');
%     PatTitle = fread(fid, 41, '*char');                 %{not used}
%     ZeroChar2 = fread(fid, 1, '*char');
%     PatComment = fread(fid, 81, '*char');               %{not used}
%     ZeroChar3 = fread(fid, 1, '*char');
%     MeasTime = fread(fid, 1, 'int32')                   %			*
%     UserId = fread(fid, 13, '*char');                   %			*
%     PatId = fread(fid, 13, '*char');                    %			*
%     Ethnic = fread(fid, 13, '*char');                   %			*
%     MeasMode = fread(fid, 13, '*char'); 				%			*
%     MacroNo = fread(fid, 13, '*char');                  %			*
%     AnaMacro = fread(fid, 13, '*char');                 %			*
%     Free = fread(fid, 174, 'int8');                     %	: Array[1..174] of Byte;					*
% 
% %PicInfoRec(Information about the ct-image)
% 
%     RecSize__	=fread(fid, 1, 'int32');                %{ = SizeOf(PicInfoRecType) }
%     PicX0 = fread(fid, 1, 'uint16')                     %{ Usage: left  edge of image}
%     PicY0 = fread(fid, 1, 'uint16')                     %{ Usage: upper edge of image}
%     PicMatrixX = fread(fid, 1, 'uint16')                %{ X-Matrix size of image    }
%     PicMatrixY = fread(fid, 1, 'uint16')                %{ Y-Matrix size of image    }
%     PicVoxelSizeX = fread(fid, 6, 'int8');              %: Real;	{ not used }
%     PicVoxelSizeY = fread(fid, 6, 'int8');              %: Real; 	{ not used }
%     Free_ = fread(fid, 64, 'int8');                     %	: Array[1..64] of Byte;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(full_name, 'r');  

fseek(fid,12,'bof');
VoxelSize = fread(fid, 1, 'double');                % Size of a voxel [mm]

fseek(fid,318,'bof');
obj_length = fread(fid, 1, 'double')                % Size of a voxel [mm]

fseek(fid,326,'bof');
percent = fread(fid, 1, 'double')                % Size of a voxel [mm]

slice_z = obj_length*percent/100

fseek(fid,1099,'bof');    
PatName = fread(fid, 41, '*char')';                 %{Patient familyname & firstname }
PatNameChar = int8(PatName(1));                     %The first byte is the number of characters used
PatName = PatName(2:(1+PatNameChar));               %Pulling the Patient Name out
    
fseek(fid,1529,'bof');
PicMatrixX = fread(fid, 1, 'uint16');               %{ X-Matrix size of image    }
PicMatrixY = fread(fid, 1, 'uint16');               %{ Y-Matrix size of image    }
    
%CT Image Data
fseek(fid,1609,'bof');
pQCT = fread(fid, PicMatrixX*PicMatrixY, 'int16');
    
%Reshape the vector to the image dimensions in the header
pQCTmat = reshape(pQCT, PicMatrixX, PicMatrixY);
   
% Close the file    
fclose(fid);

resolution = [VoxelSize VoxelSize];



A=abs(pQCTmat);                         %throw away negative air value noise
                                       



