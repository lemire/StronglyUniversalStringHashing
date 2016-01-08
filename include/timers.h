#ifndef TIMERS
#define TIMERS

typedef unsigned long long ticks;

// Taken from stackoverflow (see http://stackoverflow.com/questions/3830883/cpu-cycle-count-based-profiling-in-c-c-linux-x86-64)
// Can give nonsensical results on multi-core AMD processors.
static inline ticks rdtsc() {
    unsigned int lo, hi;
    asm volatile (
        "cpuid \n" /* serializing */
        "rdtsc"
        : "=a"(lo), "=d"(hi) /* outputs */
        : "a"(0) /* inputs */
        : "%ebx", "%ecx");
    /* clobbers*/
    return ((unsigned long long) lo) | (((unsigned long long) hi) << 32);
}

static inline ticks startRDTSC(void) {
    return rdtsc();
}

static inline ticks stopRDTSCP(void) {
    return rdtsc();
}
// start and stop are as recommended by
// Gabriele Paoloni, How to Benchmark Code Execution Times on Intel® IA-32 and IA-64 Instruction Set Architectures
// September 2010
// http://edc.intel.com/Link.aspx?id=3954

/*static __inline__ ticks fancystartRDTSC(void) {
 unsigned cycles_low, cycles_high;
 asm volatile ("CPUID\n\t"
 "RDTSC\n\t"
 "mov %%edx, %0\n\t"
 "mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::
 "%rax", "%rbx", "%rcx", "%rdx");
 return ((ticks) cycles_high << 32) | cycles_low;
 }

 static __inline__ ticks fancystopRDTSCP(void) {
 unsigned cycles_low, cycles_high;
 /// This should work fine on most machines, if the RDTSCP thing
 /// fails for you, use the  rdtsc() call instead.
 asm volatile("RDTSCP\n\t"
 "mov %%edx, %0\n\t"
 "mov %%eax, %1\n\t"
 "CPUID\n\t": "=r" (cycles_high), "=r" (cycles_low):: "%rax",
 "%rbx", "%rcx", "%rdx");
 return ((ticks) cycles_high << 32) | cycles_low;
 }*/

#endif
