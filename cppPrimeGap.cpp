// cppConsoleHelloWorld.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <immintrin.h>
//#include "InstructionSet.h"
#include <string.h>
#include <time.h>
#include <windows.h>
#include <signal.h>//#include <csignal>

#define DRILL 0
void fR(int n, int iLen, char* dst, char* primes) {
	for (int i = 0; i < iLen; i++)
	{
		dst[i] = n % primes[i];
	}
}

bool ScanForZero(int iLen, char* arr) {
	/// scan fixed length of memory adddress for null char (aka 0)
	for (int i = 0; i < iLen; i++)
	{
		if (arr[i] == 0) return false;
	}
	return true;
}
#if FALSE

inline int ScanForZero_ASM_REPSCA(int iLen, char* arr) {
	/// USES EDI,ECX,CLD
	/// input: EDI the memory address, ECX for search length, CLD
	/// output: ECX>=1 if null char is found, 0 otherwise
	__asm {
		MOV EDI, a;
		/*	LEA ESI, a ;I don't why it result ESI a different value between LEA ESI, a and  MOV EDI, a  */
		MOV ECX, iLen;
		MOV AL, 0;
		CLD;
		REPZ SCASB; //<-- core method
		JECXZ NO_FOUND;
	}
	return 1;
NO_FOUND:
	return 0;
}
bool FoundAnyDiffernt(int iLen, char* a, char* b) {
	/// compare and return true if it has some different byte between 2 memory blocks. with a fixed length.
	for (int i = 0; i < iLen; i++)
	{
		if (a[i] != b[i])
		{
			return true;
		}
	}
	return false;
}
inline int FoundAnyDiffernt_ASM_REPCMP(int iLen, char* a, char* b) {
	/// USES ESI,EDI,ECX,CLD
	/// input: EDI the memory address, ECX for search length, CLD
	/// output: ECX>=1 if null char is found, 0 otherwise
	__asm {
		MOV ESI, a;
		MOV EDI, b;
		MOV ECX, iLen;
		CLD;
		REPZ CMPSB;
		JECXZ NO_FOUND;
	}
	return 1;
NO_FOUND:
	return 0;
}
inline int FoundAnyMatch_ASM_REPCMP(int iLen, char* a, char* b) {
	__asm {
		MOV ESI, a;
		MOV EDI, b;
		MOV ECX, len;
		CLD;
		REPNE CMPSB;
		JECXZ NO_FOUND;
	}
	return 1;
NO_FOUND:
	return 0;
}
bool INC_XCARRY(int iLen, char* R, char* C) {
	for (int i = 0; i < iLen; i++)
	{
		if (R[i] < C[i])
		{
			R[i]++;
			return true;
		}
		else
		{
			R[i] = 0;
		}
	}
	return false;
}
bool DEC_XCARRY(int iLen, char* R, char* C) {
	// R-=1 in X-carry system
	// each byte as a digit,
	// each digit carry is arbitrary (<=255) rather than 10,10,10 in decimal system
	// C as carry thresholds, e.g. a prime list {2,3,5,7,...,p} or a half {1,2,4,6,... p-1}
	// return false if reach exceed minimum 0,0,0,0,
	for (int i = 0; i < iLen; i++)//for (int i = iLen - 1; i >= 0; i--)
	{
		if (R[i] == 0)
		{
			R[i] = C[i]; //reset all lower digits if in case "1000-1=0999"
		}
		else
		{
			R[i]--;
			return true;
		}
	}
	return false;
}
inline int DEC_XCARRY_ASM(int iLen, unsigned short* RC) {
	//for optimization: shuffle carry and R
	/// AL is holding R
	//  AH is holding carry
	__asm {
		//XOR EAX, EAX
		MOV ECX, iLen
		MOV ESI, RC
		MOV EDI, RC
		CLD
	}
LOOPSTART:
	__asm {
		LODSW
		DEC AL
		JS UNDERFLOW_AFTER_SUBSTRACTION
		STOSB
		JMP RET_SUCCESS
	}
UNDERFLOW_AFTER_SUBSTRACTION:
	__asm {
		MOV AL, AH
		STOSW
		LOOP LOOPSTART
	}
REACH_MIN:
	return 0;
RET_SUCCESS:
	return 1;
}
bool DEC_XCARRY_avoidZero(int iLen, char* R, char* C) {
	/// a variant of DEC_XCARRY
	/// in order to search a gap start at 1, for G2 search mission, it should always avoid 0 occurring in R sequence.
	/// return false when R={1,1,1,1...} has reached minimum limit {0,0,0,...}

	for (int i = 0; i < iLen; i++)//for (int i = iLen - 1; i >= 0; i--)
	{
		if (R[i] > 1)
		{
			R[i]--;
			return true;
		}
		else
		{
			R[i] = C[i];//reset all lower digits if in case "1000-1=0999"
		}
	}
	return false;
}

void INC_R(int iLen, char* R, char* C) {
	//increase all remainder digits by 1, this method == fR(fReverseR(R_n)+1)
	//with no carry performed
	for (int i = 0; i < iLen; i++)
	{
		if (R[i] < C[i])
		{
			R[i]++;
		}
		else
		{
			R[i] = 0;
		}
	}
}
void DEC_R(int iLen, char* R, char* C) {
	for (int i = 0; i < iLen; i++)
	{
		if (R[i] == 0)
		{
			R[i] = C[i];
		}
		else
		{
			R[i]--;
		}
	}
}
inline void DEC_R_ASM(int iLen, unsigned short* RC) {
	//for optimization: shuffle C and R
	/// AL is holding R[i]
	//  AH is holding C[i]
	__asm {
		//XOR EAX, EAX
		MOV ECX, iLen
		MOV ESI, RC
		MOV EDI, RC
		CLD
		JMP LOOPSTART
	}
LOOPSTART:
	__asm {
		LODSW
		DEC AL
		JNS LOOPTAIL
	}
UNDERFLOW_AFTER_SUBSTRACTION:
	__asm {
		MOV AL, AH
	}
LOOPTAIL:
	__asm {
		STOSB
		INC EDI
		LOOP LOOPSTART
	}
}
typedef char				I8;
typedef unsigned char		I8U;
typedef short				I16;
typedef unsigned short		I16U;
typedef int					I32;
typedef unsigned int		I32U;
typedef long long			I64;
typedef unsigned long long	I64U;

int test_DEC_XCARRY_ASM_main()
{

	//char primesdiv2plus1NoZero[] = { 1,1,3,4 };
	//char R[] = { 1,1,3,4 };
	char R[] = { 1,1,1,1,3,3,4 ,4 };
	int x;
	do
	{
		x = DEC_XCARRY_ASM(4, (unsigned short*)R);
		std::cout << R[0] << R[2] << R[4] << R[6];
	} while (x);


	char a[] = "ABCDEFG";
	//char b[] = "aBCDEFG";
	char b[] = "abcdefg";
	int len = strlen(a);
	if (FoundAnyMatch_ASM_REPCMP(7, a, b))
	{
		std::cout << "found";
	}
	else {
		std::cout << "no found";
	}
	return 0;
}

int TEST_XOR_vs_MOV_Performance() {
	//test if "XOR AX,AX" is faster than "MOV AX, 0"
	//test result: after test 100 times 10,000,000 loops
	//conclusion: no conclusive result.
	signed long long t1, t2;
	LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
	LARGE_INTEGER Frequency;

	int loopcount = 10000000;

	QueryPerformanceFrequency(&Frequency);


	QueryPerformanceCounter(&StartingTime);
	__asm MOV ECX, loopcount
	TEST1 :
	__asm MOV EAX, 0
	__asm LOOP TEST1
	QueryPerformanceCounter(&EndingTime);
	t1 = EndingTime.QuadPart - StartingTime.QuadPart;
	//printf("%llu\n", t1);

	StartingTime.QuadPart = 0;
	EndingTime.QuadPart = 0;
	QueryPerformanceCounter(&StartingTime);
	__asm MOV ECX, loopcount
	TEST2 :
	__asm XOR EAX, EAX
	__asm LOOP TEST2
	QueryPerformanceCounter(&EndingTime);
	t2 = EndingTime.QuadPart - StartingTime.QuadPart;
	//printf("%llu\n", t2);



	//
	// We now have the elapsed number of ticks, along with the
	// number of ticks-per-second. We use these values
	// to convert to the number of elapsed microseconds.
	// To guard against loss-of-precision, we convert
	// to microseconds *before* dividing by ticks-per-second.
	//
	//ElapsedMicroseconds.QuadPart = t1 - t2;
	//ElapsedMicroseconds.QuadPart *= 1000000;
	//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart
	//printf("%llu\n", ElapsedMicroseconds.QuadPart);

	return (int)(t1 - t2);
	//printf("%lld\n", t1 - t2);


}
int test_performance_main()
{
	int sum = 0, avg = 0, n = 100;
	int r;
	for (size_t i = 0; i < n; i++)
	{
		r = TEST_XOR_vs_MOV_Performance();
		printf("%d\n", r);
		sum += r;
	}
	avg = sum / n;
	printf("Average = %d", avg);
	return 0;
}

#endif // 0

static INT8 primes[] = {
2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
73, 79, 83, 89, 97, 101, 103, 107, 109,
113, 127,
//(INT8)131, 137, 139, 149, 151, 157, 163,167,
//173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
//229, 233, 239, 241, 251
};
typedef struct {
	//UINT8 p;
	//INT8 halfp;	//=p/2
	INT8 c_lb;	//the lowerbound carry of i-th element of remainder sequence in prime-nary system, =-p/2, e.g. if p=5, c_lb=-2
	INT8 c_ub;	//the upperbound carry of i-th element of remainder sequence in prime-nary system, =p/2, e.g. if p=5, c_ub=2
	INT8 n;		//i-th element of R_n, aka current position on chessboard, range [-p/2,p/2]
	INT8 absn;
	// e.g. in case of p=5, possible value are {-2,-1,0,1,2} 
	// in true form {-2,-1,0,1,2} actual represent as {128|2,128|1,0,1,2}
	INT8 rook[4];	//i-th element of R of rook_1, rook_2, .... 
	//since each element of rook_1 == -rook_2, G2_method4 only uses rook_1 whose each element in its R is absoluted.
} PACKED_ELEMENT;//packed remainder sequence of n, c and each rook config
typedef struct {
	UINT16 delta;//distance between current position n to last remaining element
	INT8 r[6]; //stores {r_1, r_2, ...} of R seq of rook config
	//for G2_method4, it only uses r_1 and r_2 as private individual part, r_3 and other elements are shared
} DIVIDED_ROOK_CONFIG;
typedef struct {
	UINT32 len_FILE;	//entire structure size<= filesize
	double timeused;
	UINT64 check_sum;//MAGIC_NUM = 0xFF0C0B0A;
	UINT32 cnt_ROOK;//main parameter of G_r(p#) function, always =2
	UINT32 maxrecord;//record the max of gap length
	UINT32 delta;//distance between current position n to last remaining element (sieved by shared subsequence of rook config. this is a optimization, so it does not need to increase gap len of each individual rook config.
	UINT32 len_RSEQ;	//main parameter of G_r(p#) function, length of remainder sequence = COUNT_PRIME_LIST-2, prime factor 2 and 3 is ignored.
	PACKED_ELEMENT SHARED[64]; //remainder sequence, each element is packed as {n, c, rook} combination
	//optimization: batch processing multiple rook configurations in one loop (which is iterating n from 0 to p#/2 once), rather than processing one rook configuration in one loop
	//e.g. in case of G2(17#)
	//there are 8 valid and distinctive 2-rook configurations (position arrangement) at dimension p=17 on chessboard (can not be gotten by shifting other arrangement and can not equal 0)
	//they are {{1,-1},{2,-2},{3,-3},{4,-4},{5,-5},{6,-6},{7,-7},{8,-8}} //-x === 17-x (mod 17)
	//if LEN_SUBSEQ=1 and BATCH_SIZE=8
	UINT32 len_DIVIDED_SUBSEQ; //length of divided subsequence, this is an optimization for batch multiple rook config
	UINT32 cnt_DIVIDED_ROOK_CONFIGS; //number of divided individual rook configs
	DIVIDED_ROOK_CONFIG DIVIDED[]; //individual rook configs
} SAVEFILE;
//int Gap_G2_sym_11(int iLen, char* R, char* R_n, char* primesMod2plus1) {
//	/// measure the gap length of certain rook position in 2-rook situation with G2sym and {r1,r2}={1,1} constraint
//	/// R as position sequence R of Both Rooks in chessboard
//	/// primesdiv2plus1NoZero
//	//***this functin is unfinished***
//	int i = 6;
//	while (FoundAnyMatch_ASM_REPCMP(iLen, R, R_n))
//	{
//		DEC_XCARRY(iLen, R, primesMod2plus1);
//		i += 6;
//		R_n += 6 * iLen;
//	};
//	return i;
//}

HANDLE hIn = GetStdHandle(STD_INPUT_HANDLE);
HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
char m_writelogbuffer[512];

int cnt_THREAD;
int len_RSEQ;
int cnt_DIVIDED_ROOK_CONFIGS;
int len_FILE;
int len_DIVIDED_SUBSEQ;
PACKED_ELEMENT* R;
int idx_RSEQ_FIRST = 0;
int idx_RSEQ_LAST_DIVIDED;
int idx_RSEQ_FIRST_SHARED;
int idx_RSEQ_LAST;
int idx_DIVIDED_LAST;
CRITICAL_SECTION cs;
SYNCHRONIZATION_BARRIER barrier;
volatile LONG m_loopcounter;

char lpFilenameSave[] = "G2_method4.save";
HANDLE m_hFileSave;
char lpFilenameLog[] = "G2_method4(###).log";
HANDLE m_hFileLog;
SAVEFILE* m_lpSaveFile;
char m_random_text_buffer[512];

time_t timeStart;
time_t timeEnd;
enum WORK_MESSAGE
{
	DISMISS,
	NOP,
	PROCESS_MATCHED_ON_END_LOOP_N,
	PROCESS_NOT_MATCHED_ON_END_LOOP_N,
	PROCESS_NOT_MATCHED_DURING_LOOP_N,
};
volatile WORK_MESSAGE wMsg;

void WriteLog(char* text)
{
	char* p = m_writelogbuffer;
	size_t l = sizeof(m_writelogbuffer); //len_remaining
	size_t i;
	size_t len = 0; //nNumberOfBytesToWrite

	time_t rawtime;
	time(&rawtime);
	struct tm info;
	gmtime_s(&info, &rawtime);
	i = strftime(p, l, "%Y-%m-%d %H:%M:%SZ\t", &info);p += i;l -= i;len += i;

	i = sprintf_s(p, l, "%s\n", text);p += i;l -= i;len += i;
	//*p = '\n';//append '\n'

	{
		DWORD NumberOfBytesWritten;
		WriteConsoleA(hOut,
			m_writelogbuffer,
			(DWORD)len,
			&NumberOfBytesWritten,
			NULL);
#if DRILL
#else
		WriteFile(
			m_hFileLog,
			m_writelogbuffer,
			(DWORD)len,
			&NumberOfBytesWritten,
			NULL
		);
#endif
	}
}

void LogNewRecord(DIVIDED_ROOK_CONFIG* rookconfig, bool isNewRecord, bool isCrossedHalfChessboard) {
	char* p = m_writelogbuffer;
	size_t l = sizeof(m_writelogbuffer); //len_remaining
	size_t i;
	size_t len = 0; //nNumberOfBytesToWrite

	int j;

	time_t rawtime;
	time(&rawtime);
	struct tm info;
	gmtime_s(&info, &rawtime);
	i = strftime(p, l, "%Y-%m-%d %H:%M:%SZ\t", &info);p += i;l -= i;len += i;
	i = sprintf_s(p, l, "{gaplen:%d,R_rook:{", m_lpSaveFile->maxrecord);p += i;l -= i;len += i;//prefix {1,1} for prime factor 2 and 3

	for (j = idx_RSEQ_LAST; j >= idx_RSEQ_FIRST_SHARED; j--) { //5,7,11 .. print shared part
		i = sprintf_s(p, l, "%d,", *R[j].rook);p += i;l -= i;len += i;
	}
	for (j = idx_RSEQ_LAST_DIVIDED; j >= 0; j--) {//print unshared part
		i = sprintf_s(p, l, "%d,", *(rookconfig->r + j));p += i;l -= i;len += i;
	}

	i = sprintf_s(p, l, "},R_n:{");p += i;l -= i;len += i;
	for (j = idx_RSEQ_LAST; j >= 0; j--)
	{
		i = sprintf_s(p, l, "%d,", R[j].n);p += i;l -= i;len += i;
	}
	i = sprintf_s(p, l, "}}");p += i;l -= i;len += i;

	if (isNewRecord)
	{
		i = sprintf_s(p, l, ", new record");p += i;l -= i;len += i;
	}
	if (isCrossedHalfChessboard)
	{
		i = sprintf_s(p, l, ", cross p#/2");p += i;l -= i;len += i;
	}
	i = sprintf_s(p, l, "\n");p += i;l -= i;len += i;
	{
		DWORD NumberOfBytesWritten;
		WriteConsoleA(hOut,
			m_writelogbuffer,
			(DWORD)len,
			&NumberOfBytesWritten,
			NULL);
#if DRILL
#else
		WriteFile(
			m_hFileLog,
			m_writelogbuffer,
			(DWORD)len,
			&NumberOfBytesWritten,
			NULL
		);
#endif // DRILL
	}
}

void DebugOutput() {
	char* p = m_writelogbuffer;
	size_t l = sizeof(m_writelogbuffer); //len_remaining
	size_t i;
	size_t len = 0; //nNumberOfBytesToWrite
	int j;

	//i = sprintf_s(p, l, "R_n={");p += i;l -= i;len += i;
	//for (j = idx_RSEQ_LAST; j >= 0; j--)
	//{
	//	i = sprintf_s(p, l, "%d,", R[j].n);p += i;l -= i;len += i;
	//}
	//i = sprintf_s(p, l, "}");p += i;l -= i;len += i;

	//i = sprintf_s(p, l, "delta={");p += i;l -= i;len += i;
	//for (int idxConfig = 0; idxConfig <= idx_DIVIDED_LAST; idxConfig++)
	//{
	//	DIVIDED_ROOK_CONFIG* rookconfig = m_lpSaveFile->DIVIDED + idxConfig;
	//	i = sprintf_s(p, l, "%d,", rookconfig->delta);p += i;l -= i;len += i;
	//}
	//i = sprintf_s(p, l, "}\n");p += i;l -= i;len += i;

	i = sprintf_s(p, l, "{maxrecord:%d,R_rook:{", m_lpSaveFile->maxrecord);p += i;l -= i;len += i;//prefix {1,1} for prime factor 2 and 3

	for (j = idx_RSEQ_LAST; j >= idx_RSEQ_FIRST_SHARED; j--) { //5,7,11 .. print shared part
		i = sprintf_s(p, l, "%d,", *R[j].rook);p += i;l -= i;len += i;
	}
	for (j = idx_RSEQ_LAST_DIVIDED; j >= 0; j--) {//print unshared part
		i = sprintf_s(p, l, "X,");p += i;l -= i;len += i;
	}
	i = sprintf_s(p, l, "}}\n");p += i;l -= i;len += i;


	DWORD NumberOfBytesWritten;
	WriteConsoleA(hOut,
		m_writelogbuffer,
		(DWORD)len,
		&NumberOfBytesWritten,
		NULL);
}
inline void TryUpdateGlobalRecord(DIVIDED_ROOK_CONFIG* rookconfig, bool isCrossedHalfChessboard) {
	EnterCriticalSection(&cs);
	if (rookconfig->delta < m_lpSaveFile->maxrecord) {

	}
	else if (rookconfig->delta > m_lpSaveFile->maxrecord) {
		m_lpSaveFile->maxrecord = (UINT32)rookconfig->delta;
		LogNewRecord(rookconfig, true, isCrossedHalfChessboard);
	}
	else //if (rookconfig->delta == m_lpSaveFile->maxrecord) {
	{
		LogNewRecord(rookconfig, false, isCrossedHalfChessboard);
	}
	LeaveCriticalSection(&cs);
}
DWORD WINAPI THREAD_START_ROUTINE_G2_method4_WorkerThread(__in LPVOID lpParameter) {
	const int id_THREAD = (int)lpParameter;
	int idxConfig;
	int idxSeq;
	do
	{
		EnterSynchronizationBarrier(&barrier, 0); // SYNCHRONIZATION_BARRIER_FLAGS_NO_DELETE
		switch (wMsg)
		{
		default:
		case NOP:
			EnterCriticalSection(&cs);
			printf("ThreadId: %d NOPed\n", id_THREAD);
			LeaveCriticalSection(&cs);
			break;
		case DISMISS:
			EnterCriticalSection(&cs);
			printf("ThreadId: %d exit\n", id_THREAD);
			LeaveCriticalSection(&cs);
			return 0;//ExitThread(0);
			break;
		case PROCESS_MATCHED_ON_END_LOOP_N:
			for (idxConfig = id_THREAD; idxConfig <= idx_DIVIDED_LAST; idxConfig += cnt_THREAD)
			{
				DIVIDED_ROOK_CONFIG* rookconfig = m_lpSaveFile->DIVIDED + idxConfig;
				rookconfig->delta += m_lpSaveFile->delta;
				rookconfig->delta *= 2;
				rookconfig->delta--;
				TryUpdateGlobalRecord(rookconfig, true);
			}
			break;
		case PROCESS_NOT_MATCHED_ON_END_LOOP_N:
			for (idxConfig = id_THREAD; idxConfig <= idx_DIVIDED_LAST; idxConfig += cnt_THREAD)
			{
				DIVIDED_ROOK_CONFIG* rookconfig = m_lpSaveFile->DIVIDED + idxConfig;
				INT8* S = rookconfig->r;
				for (idxSeq = idx_RSEQ_LAST_DIVIDED; idxSeq >= 0; idxSeq--)
				{
					if (S[idxSeq] == R[idxSeq].absn)
					{
						//now n=p#/2+1, if not matched, n=p#/2 also not match
						rookconfig->delta *= 2;
						rookconfig->delta++;
						TryUpdateGlobalRecord(rookconfig, true);
						break;
					}
				}
			}
			break;
		case PROCESS_NOT_MATCHED_DURING_LOOP_N:
			for (idxConfig = id_THREAD; idxConfig <= idx_DIVIDED_LAST; idxConfig += cnt_THREAD)
			{
				DIVIDED_ROOK_CONFIG* rookconfig = m_lpSaveFile->DIVIDED + idxConfig;
				INT8* S = rookconfig->r;
				rookconfig->delta += m_lpSaveFile->delta;
				for (idxSeq = idx_RSEQ_LAST_DIVIDED; idxSeq >= 0; idxSeq--)
				{
					if (S[idxSeq] == R[idxSeq].absn)
					{
						//rookconfig->delta += m_lpSaveFile->delta;
						goto NEXT_CONFIG;
					}
				}
				//rookconfig->delta += m_lpSaveFile->delta;
				TryUpdateGlobalRecord(rookconfig, false);
				rookconfig->delta = 0;
			NEXT_CONFIG:;
			}
			break;
		}
		//if (InterlockedDecrement(&m_cntWorkersAtWork) == 0) //Sometime InterlockedDecrement misscounts 1 causing entire program stall. it is not suitable for thread sync
		//{
		//	//InterlockedCompareExchange(&m_cntWorkersAtWork, cnt_THREAD, 1);
		//	InterlockedExchange(&m_cntWorkersAtWork, cnt_THREAD);
		//	ResetEvent(m_hEventStartWork);
		//	SetEvent(m_hEventWorkDone);
		//}
		//else
		//{
		//	WaitForSingleObject(m_hEventWorkDone, INFINITE);
		//}
		EnterSynchronizationBarrier(&barrier, 0);
	} while (id_THREAD);
}
void WorkersDoWork(WORK_MESSAGE msg) {
	wMsg = msg;
	THREAD_START_ROUTINE_G2_method4_WorkerThread(0);
}
UINT64 CheckSum64(size_t filesize, void* lpbaseaddress) {
	UINT64* ptr = (UINT64*)lpbaseaddress;
	UINT64* ptrEnd = (UINT64*)((char*)lpbaseaddress + filesize);
	UINT64 sum = 0;
	for (; ptr < ptrEnd; ptr++)
	{
		sum += *ptr;
	}
	return sum;
}
void SaveSaveFile() {
	{
		m_lpSaveFile->check_sum = 0;
		UINT64 sum = CheckSum64(len_FILE, m_lpSaveFile);
		UINT64 ComplementSum = 0 - sum;
		m_lpSaveFile->check_sum = ComplementSum;
	}
	//ROLLBACK_MAPVIEW:
	//bFlag = 
	FlushViewOfFile(m_lpSaveFile, 0);
}
void G2_method4(bool isNewStart, bool isUsingFileMapping) {
	R = m_lpSaveFile->SHARED;
	idx_RSEQ_FIRST = 0;
	idx_RSEQ_LAST_DIVIDED = m_lpSaveFile->len_DIVIDED_SUBSEQ - 1;
	idx_RSEQ_FIRST_SHARED = m_lpSaveFile->len_DIVIDED_SUBSEQ;
	idx_RSEQ_LAST = m_lpSaveFile->len_RSEQ - 1;
	idx_DIVIDED_LAST = m_lpSaveFile->cnt_DIVIDED_ROOK_CONFIGS - 1;
	if (cnt_THREAD > cnt_DIVIDED_ROOK_CONFIGS)
	{
		cnt_THREAD = cnt_DIVIDED_ROOK_CONFIGS;
	}
	InitializeCriticalSection(&cs);
	InitializeSynchronizationBarrier(&barrier, cnt_THREAD, -1);
	printf("Creating %d thread(s)...\n", cnt_THREAD);
	HANDLE* hThreadArray = (HANDLE*)malloc(sizeof(HANDLE) * cnt_THREAD);
	for (size_t i = 1; i < cnt_THREAD; i++)
	{
		hThreadArray[i] = CreateThread(NULL,
			0,
			THREAD_START_ROUTINE_G2_method4_WorkerThread,
			(LPVOID*)i,
			0,
			NULL);
		if (hThreadArray[i] == NULL)
		{
			printf("Error while creating thread\n");
			return;
		}
	}
	Sleep(1000);

	int i;
	int idxConfig;
	int idxSeq;
	INT8* n, * absn, * r;
	COORD coordZero;
	coordZero.X = 0; coordZero.Y = 0;
	CONSOLE_SCREEN_BUFFER_INFO bufferInfo;
	COORD* coordCurrent = &bufferInfo.dwCursorPosition;

	if (isNewStart)
	{
		//INITIALIZE_RSEQ_SAVEFILE:
		for (i = idx_RSEQ_LAST; i >= 0; i--) //invert order {...,17,13,11,7,5}
		{
			int j = idx_RSEQ_LAST - i + 2;//Rseq[0].p = 5, so i=0, prime[2]
			INT8 p = primes[j];
			INT8 half_p = p / 2;
			R[i].c_lb = -half_p;
			R[i].c_ub = half_p;
			R[i].n = 0;
			R[i].absn = 0;
			*R[i].rook = half_p;
		}
		DIVIDED_ROOK_CONFIG* rookconfig = m_lpSaveFile->DIVIDED + idx_DIVIDED_LAST;
		INT8* S = rookconfig->r;
		for (idxSeq = idx_RSEQ_LAST_DIVIDED; idxSeq >= 0; idxSeq--)
		{
			S[idxSeq] = R[idxSeq].c_ub;
		}
		for (idxConfig = idx_DIVIDED_LAST - 1; idxConfig >= 0; idxConfig--)
		{
			*(rookconfig - 1) = *rookconfig;
			rookconfig--;
			S = rookconfig->r;
			for (idxSeq = 0; idxSeq <= idx_RSEQ_LAST_DIVIDED; idxSeq++) //DEC_XCARRY_FOR_ROOKCONFIG
			{
				S[idxSeq]--;
				if (S[idxSeq] == 0)
				{
					S[idxSeq] = R[idxSeq].c_ub; //reset digits in case "1000-1=0999"
				}
				else
				{
					break;
				}
			}
		}
	}

	LONG loopcounter = 0;
	//LONG loop_n = 1;
	//LONG loop_r = 1;
	//for (i = idx_RSEQ_LAST; i >= 0; i--) //loop count n
	//{
	//	int j = idx_RSEQ_LAST - i + 2;
	//	INT8 p = primes[j];

	//	loop_n *= p;
	//}
	//loop_n /= 2;
	//loop_n++;
	//for (i = idx_RSEQ_LAST; i >= idx_RSEQ_FIRST_SHARED; i--)//loop count r
	//{
	//	loop_r *= R[i].c_ub;
	//}
	//loopcounter = loop_n * loop_r;
	//loopcounter++;
	//loopcounter *= loopcounter;
	//InterlockedExchange(&m_loopcounter, loopcounter);
	InterlockedExchange(&m_loopcounter, 0);
	{
	NEXT_N:
		m_lpSaveFile->delta++;

		for (i = idx_RSEQ_LAST; i >= 0; i--)//n++ until n reaches p#/2
		{
			n = &R[i].n;
			absn = &R[i].absn;
			if (*n == R[i].c_ub)
			{
				*n = R[i].c_lb;
			}
			else
			{
				(*n)++;
				*absn = *n < 0 ? -*n : *n;
				goto N_LOOP_NOT_END;
			}
		}

		//{
		//	//SetConsoleCursorPosition(hOut, coord);
		//	//printf("n=%d,", temp);
		//	DebugOutput();
		//}
	//N_LOOP_END: //n == p#/2, has tried every value of n for for G2_method4, one loop ended
		for (i = idx_RSEQ_LAST; i >= idx_RSEQ_FIRST_SHARED; i--)
		{
			if (R[i].absn == *R[i].rook)
			{
				WorkersDoWork(WORK_MESSAGE::PROCESS_MATCHED_ON_END_LOOP_N);
				goto NEXT_R;
			}
		}
		//N_LOOP_END_AND_NO_MATCH_IN_SHARED_PART:

		WorkersDoWork(WORK_MESSAGE::PROCESS_NOT_MATCHED_ON_END_LOOP_N);
		goto NEXT_R;

	N_LOOP_NOT_END:
		for (i--; i >= 0; i--)//continue increase rest element of R_n
		{
			n = &R[i].n;
			absn = &R[i].absn;
			if (*n == R[i].c_ub)
			{
				*n = R[i].c_lb;
			}
			else
			{
				(*n)++;
				*absn = *n < 0 ? -*n : *n;
			}
		}

		//{
		//	//SetConsoleCursorPosition(hOut, coord);
		//	//printf("n=%d,", temp);
		//	DebugOutput();
		//}

		for (i = idx_RSEQ_LAST; i >= idx_RSEQ_FIRST_SHARED; i--) //check match if rook=n
		{
			if (R[i].absn == *R[i].rook)
			{
				//N_LOOP_NOT_END_AND_MATCH_IN_SHARED_PART:
				goto NEXT_N;
			}
		}
		//N_LOOP_NOT_END_AND_NO_MATCH_IN_SHARED_PART:
		WorkersDoWork(WORK_MESSAGE::PROCESS_NOT_MATCHED_DURING_LOOP_N);
		m_lpSaveFile->delta = 0;
		goto NEXT_N;
	}
	{
	NEXT_R://iterate and try next rook config (decrease shared part)
		m_lpSaveFile->delta = 0; //reset		
		for (i = idx_RSEQ_LAST; i >= 0; i--)//reset n
		{
			R[i].absn = R[i].n = 0;
		}
		for (idxConfig = 0; idxConfig <= idx_DIVIDED_LAST; idxConfig++)
		{
			DIVIDED_ROOK_CONFIG* rookconfig = m_lpSaveFile->DIVIDED + idxConfig;
			rookconfig->delta = 0;
		}
		for (i = idx_RSEQ_FIRST_SHARED; i <= idx_RSEQ_LAST; i++)
		{
			r = R[i].rook;
			(*r)--;
			if (*r == 0)
			{
				*r = R[i].c_ub; //reset digits in case "1000-1=0999"
			}
			else
			{
				////loop R not end //continue increase rest element of R_rook
				//for (i++; i <= idx_RSEQ_LAST; i++)
				//{
				//	r = R[i].rook;
				//	(*r)--;
				//	if (*r == 0)
				//	{
				//		*r = R[i].c_ub; //reset digits in case "1000-1=0999"
				//	}
				//}
				{
					GetConsoleScreenBufferInfo(hOut, &bufferInfo);
					SetConsoleCursorPosition(hOut, *coordCurrent);
					DebugOutput();
					SetConsoleCursorPosition(hOut, *coordCurrent);
				}
#if DRILL
#else
#endif // DRILL
				//loopcounter = InterlockedDecrement(&m_loopcounter);
				loopcounter = InterlockedIncrement(&m_loopcounter);
				if (loopcounter == 0)
				{
					//program is interrupted, may be Ctrl+C
					WriteLog((char*)"paused.");
					goto FINALIZE_WORKER_THREADS;
				}
				//{
				//	GetConsoleScreenBufferInfo(hOut, &bufferInfo);
				//	SetConsoleCursorPosition(hOut, *coordCurrent);
				//	printf("[ %0000x ]", loopcounter);
				//	SetConsoleCursorPosition(hOut, *coordCurrent);
				//}

				//if (isUsingFileMapping)
				//{
				//	SaveSaveFile();
				//}
				goto NEXT_N;
			}
		}
		//now rook config == {1,1,1 ...}, has tried every possible R of rook config, search ends
		goto FINALIZE_WORKER_THREADS;
	}
FINALIZE_WORKER_THREADS:
	WorkersDoWork(DISMISS);//DismissWorkers();
	WaitForMultipleObjects(cnt_THREAD - 1, hThreadArray + 1, TRUE, INFINITE);
	for (int i = 1; i < cnt_THREAD; i++)
	{
		CloseHandle(hThreadArray[i]);
	}
	free(hThreadArray);
	DeleteSynchronizationBarrier(&barrier);
	DeleteCriticalSection(&cs);
	return;
}
bool isNewStart = true;
bool isUsingFileMapping = false;
BOOL bFlag;   // a result holder
HANDLE hFileMappingObject;      // handle for the file's memory-mapped region
DWORD dwBytesWritten;  // number of bytes written
DWORD dwFileSize = 0;//dwSysGran > 4096 ? dwSysGran : 4096;     // temporary storage for file sizes
DWORD dwFileMapSize = 0;  // size of the file mapping
DWORD dwMapViewFileOffset = 0; // where to start the file map view
DWORD dwMapViewSize = 0;  // the size of the view
LPVOID lpMapAddress;  // pointer to the base address of the
int G2_method4_main() {
	int ret = EXIT_FAILURE;
	SYSTEM_INFO SysInfo;  // system information; used to get granularity
	GetSystemInfo(&SysInfo);
	const DWORD dwSysGran = SysInfo.dwAllocationGranularity;      // system allocation granularity
	cnt_THREAD = SysInfo.dwNumberOfProcessors;


	int i;
	char c;

	{
		len_DIVIDED_SUBSEQ = 2;
		cnt_DIVIDED_ROOK_CONFIGS = 1;
		int d;
		printf("Prime Gap Search Program For My Math Paper\n");
		printf(" -- function G2_method4(pi)\n");
		printf(" by n43e120 2019/12\n");
		printf("Pleas input param 'pi' (the count of prime factors, between 3 and 10)?: ");
		scanf_s("%d", &d);
		if (d > 2 && d <= 10)
		{
			len_RSEQ = d - 2;
			if (len_RSEQ < len_DIVIDED_SUBSEQ)
			{
				len_DIVIDED_SUBSEQ = len_RSEQ;
			}
		}
		else
		{
			printf("Input out of range");
			goto ROLLBACK_0;
		}
		for (i = len_DIVIDED_SUBSEQ - 1; i >= 0; i--)//use last 2 prime factors for divided part of rook config
		{
			cnt_DIVIDED_ROOK_CONFIGS *= primes[len_RSEQ + 1 - i] / 2;
		}
		len_FILE = sizeof(SAVEFILE) + sizeof(DIVIDED_ROOK_CONFIG) * cnt_DIVIDED_ROOK_CONFIGS;
		if ((len_FILE & 63) > 0)
		{
			len_FILE = (len_FILE & ~63) + 64;
		}
	}
	{
		sprintf_s(lpFilenameLog, sizeof(lpFilenameLog), "G2_method4(%d).log", len_RSEQ + 2);
	}

	m_hFileLog = CreateFileA(lpFilenameLog,
		FILE_APPEND_DATA,
		FILE_SHARE_READ,
		NULL,
		OPEN_ALWAYS,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	if (m_hFileLog == INVALID_HANDLE_VALUE)
	{
		m_hFileLog = NULL;
		printf("Can not open file: %s\n", lpFilenameLog);
		//MessageBoxA(NULL, "Can not open file!", "Error", MB_OK);
		goto ROLLBACK_0;
	}

	printf("Run it without saving? [Y/n]");
	scanf_s(" %c", &c, 1);
	if ((('y' - 'Y') | c) == 'y')
	{
		goto NEW_START_USING_MEMORY;
	}

	isUsingFileMapping = true;
	m_hFileSave = CreateFileA(lpFilenameSave,
		GENERIC_WRITE | GENERIC_READ,
		FILE_SHARE_READ,
		NULL,
		OPEN_EXISTING,
		NULL,
		NULL);
	if (m_hFileSave == INVALID_HANDLE_VALUE) // file not exist
	{
		//printf("File: %s not exists, run it in memory anyway? [Y/n]", lpFilenameSave);
		//scanf_s(" %c", &c, 1);
		//if ((('y' - 'Y') | c) == 'y')
		//{
		//	goto NEW_START_USING_MEMORY;
		//}
		m_hFileSave = CreateFileA(lpFilenameSave,
			GENERIC_WRITE | GENERIC_READ,
			FILE_SHARE_READ,
			NULL,
			OPEN_ALWAYS,
			NULL,
			NULL);
		if (m_hFileSave == INVALID_HANDLE_VALUE)
		{
			m_hFileSave = NULL;
			printf("Can not open file: %s\n", lpFilenameSave);
			goto ROLLBACK_LOGFILE;
		}
		//goto FILE_EXISTED;
	}
	//FILE_EXISTED:
	{

		DWORD filesize = GetFileSize(m_hFileSave, NULL);
		dwFileMapSize = (DWORD)len_FILE;
		dwMapViewSize = (DWORD)len_FILE;
		hFileMappingObject = CreateFileMapping(m_hFileSave,          // current file handle
			NULL,           // default security
			PAGE_READWRITE, // read/write permission
			0,              // size of mapping object, high
			dwFileMapSize,  // size of mapping object, low
			NULL);          // name of mapping object
		if (hFileMappingObject == NULL)
		{
			printf("hMapFile is NULL: last error: %d\n", GetLastError());
			goto ROLLBACK_SAVEFILE;
		}

		// Map the view and test the results.
		lpMapAddress = MapViewOfFile(hFileMappingObject,
			FILE_MAP_ALL_ACCESS,
			0,
			dwMapViewFileOffset,
			dwMapViewSize);
		if (lpMapAddress == NULL)
		{
			printf("lpMapAddress is NULL: last error: %d\n", GetLastError());
			goto ROLLBACK_MAPFILE;
		}
		m_lpSaveFile = (SAVEFILE*)lpMapAddress;
		if (len_FILE <= filesize)
		{
			if (m_lpSaveFile->len_FILE == len_FILE)
			{
				if (CheckSum64(m_lpSaveFile->len_FILE, m_lpSaveFile) == 0)
				{
					isNewStart = false;
					goto RESUME_A_SAVE;
				} else{
					printf("savefile checksum failed. It may exit before save checksum at last run.\n");
					printf("Use it to Resume?[Y] or Exit Program?[N]:");
					scanf_s(" %c", &c, 1);
					if ((('y' - 'Y') | c) == 'y')
					{
						isNewStart = false;
						goto RESUME_A_SAVE;
					}
				}
			}
		}
		printf("savefile format is corrupted or wrong, be left for another search mission.\n");
		printf("Overwrite to a New Start?[Y] or Exit Program?[N]:");
		scanf_s(" %c", &c, 1);
		if ((('y' - 'Y') | c) == 'y')
		{
			isNewStart = true;
		}
		else
		{
			goto ROLLBACK_MAPVIEW_NO_FLUSH;
		}
	}

	//initialize task package
	if (isNewStart)//new empty file
	{
		if (!isUsingFileMapping)
		{
		NEW_START_USING_MEMORY:
			//	goto USING_MEMORY;
			//USING_MEMORY:
			m_lpSaveFile = (SAVEFILE*)malloc(len_FILE);
			if (m_lpSaveFile == NULL)
			{
				printf("malloc failed.");
				goto ROLLBACK_0;
			}
		}
		ZeroMemory(m_lpSaveFile, len_FILE);
		m_lpSaveFile->len_FILE = len_FILE;
		m_lpSaveFile->cnt_ROOK = 2;
		m_lpSaveFile->len_RSEQ = len_RSEQ;
		m_lpSaveFile->len_DIVIDED_SUBSEQ = len_DIVIDED_SUBSEQ;
		m_lpSaveFile->cnt_DIVIDED_ROOK_CONFIGS = cnt_DIVIDED_ROOK_CONFIGS;
		m_lpSaveFile->maxrecord = len_RSEQ * 2 + 1;

		//c = getchar();
#if DRILL
		c = 'y';
		//printf("\n");
#else
		printf("Start New, Ready? [Y/N]");
		scanf_s(" %c", &c, 1);
#endif // DRILL

		if ((('y' - 'Y') | c) == 'y')
		{
			{
				sprintf_s(m_random_text_buffer, sizeof(m_random_text_buffer), "Search started G2_method4(%d).", m_lpSaveFile->len_RSEQ + 2);
				WriteLog(m_random_text_buffer);
			}
			goto BEGIN_MAIN_FUNC;
		}
	}
	else
	{
	RESUME_A_SAVE:
		printf("Resume, Ready? [Y/N]");
		Sleep(1000);
		char c;
		//c = getchar();
		scanf_s(" %c", &c, 1);
		if ((('y' - 'Y') | c) == 'y')
		{
			{
				sprintf_s(m_random_text_buffer, sizeof(m_random_text_buffer), "resumed.");
				WriteLog(m_random_text_buffer);
			}
			goto BEGIN_MAIN_FUNC;
		}
	}
	goto FINALIZE_FILE_HANDLES;
	{
	BEGIN_MAIN_FUNC:
		time(&timeStart);
		G2_method4(isNewStart, isUsingFileMapping);
		time(&timeEnd);
		m_lpSaveFile->timeused += difftime(timeEnd, timeStart);
		if (InterlockedOr(&m_loopcounter, 0) != 0)
		{
			//SEARCH_SUCCESSFULLY_END:
			sprintf_s(m_random_text_buffer, sizeof(m_random_text_buffer), "Search ended, G2(%d)/6=%d, time cost=%.1f seconds", m_lpSaveFile->len_RSEQ + 2, m_lpSaveFile->maxrecord, m_lpSaveFile->timeused);
			WriteLog(m_random_text_buffer);
			ret = EXIT_SUCCESS;
		}
		else
		{
			printf("Ctrl+C pressed\n");
		}
	}
FINALIZE_FILE_HANDLES:
	if (isUsingFileMapping)
	{
		SaveSaveFile();

	ROLLBACK_MAPVIEW_NO_FLUSH:
		bFlag = UnmapViewOfFile(lpMapAddress);
	ROLLBACK_MAPFILE:
		bFlag = CloseHandle(hFileMappingObject); // close the file mapping object
		if (!bFlag)
		{
			printf("Error %ld occurred closing the mapping object!\n", GetLastError());
		}
	ROLLBACK_SAVEFILE:
		bFlag = CloseHandle(m_hFileSave);
		if (!bFlag)
		{
			printf("Error %ld occurred closing the file!\n", GetLastError());
		}
	}
	else
	{
		free((void*)m_lpSaveFile);
	}

ROLLBACK_LOGFILE:
	bFlag = CloseHandle(m_hFileLog);
	if (!bFlag)
	{
		printf("Error %ld occurred closing the file!\n", GetLastError());
	}
ROLLBACK_0:
	printf("Program exited with code %d\n.", ret);
	return ret;
}
void Test_CPU_Support() {
	//auto& outstream = std::cout;

//auto support_message = [&outstream](std::string isa_feature, bool is_supported) {
//	outstream << isa_feature << (is_supported ? " supported" : " not supported") << std::endl;
//};
//bool bSupport = InstructionSet::AVX2();
//support_message("AVX2", bSupport);
//if (!bSupport)
//{
//	return 0;
//}
//bSupport = InstructionSet::AVX512CD();
//support_message("AVX512CD", bSupport);
//if (!bSupport)
//{
//	return 0;
//}
}
BOOL FileExists(LPCTSTR szPath)
{
	DWORD dwAttrib = GetFileAttributes(szPath);

	return (dwAttrib != INVALID_FILE_ATTRIBUTES &&
		!(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
}
BOOL FileExistsA(LPCSTR szPath)
{
	DWORD dwAttrib = GetFileAttributesA(szPath);

	return (dwAttrib != INVALID_FILE_ATTRIBUTES &&
		!(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
}
COORD GetCursorPosition() {
	HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO bufferInfo;
	GetConsoleScreenBufferInfo(h, &bufferInfo);
	return bufferInfo.dwCursorPosition;
}
void SetCursorPosition(int XPos, int YPos) {
	COORD coord;
	coord.X = XPos; coord.Y = YPos;
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord);
}
void sighandler(int signum)
{
	InterlockedExchange(&m_loopcounter, -1);
	printf("Signal %d\n", signum);
}
int main()
{
	signal(SIGINT, sighandler);
	G2_method4_main();
	return 0;
}