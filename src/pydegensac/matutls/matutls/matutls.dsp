# Microsoft Developer Studio Project File - Name="matutls" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=matutls - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "matutls.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "matutls.mak" CFG="matutls - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "matutls - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "matutls - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "matutls - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "../../../lib/win32/Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "matutls - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "../../../lib/win32/Debug/"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "matutls - Win32 Release"
# Name "matutls - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\atou1.c
# End Source File
# Begin Source File

SOURCE=..\atovm.c
# End Source File
# Begin Source File

SOURCE=..\chouse.c
# End Source File
# Begin Source File

SOURCE=..\chousv.c
# End Source File
# Begin Source File

SOURCE=..\cmattr.c
# End Source File
# Begin Source File

SOURCE=..\cmcpy.c
# End Source File
# Begin Source File

SOURCE=..\cminv.c
# End Source File
# Begin Source File

SOURCE=..\cmmul.c
# End Source File
# Begin Source File

SOURCE=..\cmmult.c
# End Source File
# Begin Source File

SOURCE=..\cmprt.c
# End Source File
# Begin Source File

SOURCE=..\csolv.c
# End Source File
# Begin Source File

SOURCE=..\cvmul.c
# End Source File
# Begin Source File

SOURCE=..\eigen.c
# End Source File
# Begin Source File

SOURCE=..\eigval.c
# End Source File
# Begin Source File

SOURCE=..\evmax.c
# End Source File
# Begin Source File

SOURCE=..\hconj.c
# End Source File
# Begin Source File

SOURCE=..\heigval.c
# End Source File
# Begin Source File

SOURCE=..\heigvec.c
# End Source File
# Begin Source File

SOURCE=..\hevmax.c
# End Source File
# Begin Source File

SOURCE=..\hmgen.c
# End Source File
# Begin Source File

SOURCE=..\house.c
# End Source File
# Begin Source File

SOURCE=..\housev.c
# End Source File
# Begin Source File

SOURCE=..\ldumat.c
# End Source File
# Begin Source File

SOURCE=..\ldvmat.c
# End Source File
# Begin Source File

SOURCE=..\matprt.c
# End Source File
# Begin Source File

SOURCE=..\mattr.c
# End Source File
# Begin Source File

SOURCE=..\mcopy.c
# End Source File
# Begin Source File

SOURCE=..\minv.c
# End Source File
# Begin Source File

SOURCE=..\mmul.c
# End Source File
# Begin Source File

SOURCE=..\ortho.c
# End Source File
# Begin Source File

SOURCE=..\otrma.c
# End Source File
# Begin Source File

SOURCE=..\otrsm.c
# End Source File
# Begin Source File

SOURCE=..\psinv.c
# End Source File
# Begin Source File

SOURCE=..\qrbdi.c
# End Source File
# Begin Source File

SOURCE=..\qrbdu1.c
# End Source File
# Begin Source File

SOURCE=..\qrbdv.c
# End Source File
# Begin Source File

SOURCE=..\qrecvc.c
# End Source File
# Begin Source File

SOURCE=..\qreval.c
# End Source File
# Begin Source File

SOURCE=..\qrevec.c
# End Source File
# Begin Source File

SOURCE=..\qrlsq.c
# End Source File
# Begin Source File

SOURCE=..\rmmult.c
# End Source File
# Begin Source File

SOURCE=..\ruinv.c
# End Source File
# Begin Source File

SOURCE=..\smgen.c
# End Source File
# Begin Source File

SOURCE=..\solvps.c
# End Source File
# Begin Source File

SOURCE=..\solvru.c
# End Source File
# Begin Source File

SOURCE=..\solvtd.c
# End Source File
# Begin Source File

SOURCE=..\sv2u1v.c
# End Source File
# Begin Source File

SOURCE=..\sv2uv.c
# End Source File
# Begin Source File

SOURCE=..\sv2val.c
# End Source File
# Begin Source File

SOURCE=..\svdu1v.c
# End Source File
# Begin Source File

SOURCE=..\svduv.c
# End Source File
# Begin Source File

SOURCE=..\svdval.c
# End Source File
# Begin Source File

SOURCE=..\trncm.c
# End Source File
# Begin Source File

SOURCE=..\trnm.c
# End Source File
# Begin Source File

SOURCE=..\unitary.c
# End Source File
# Begin Source File

SOURCE=..\utrncm.c
# End Source File
# Begin Source File

SOURCE=..\utrnhm.c
# End Source File
# Begin Source File

SOURCE=..\vmul.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\complex.h
# End Source File
# Begin Source File

SOURCE=..\matutl.h
# End Source File
# End Group
# End Target
# End Project
