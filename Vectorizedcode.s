#define STDOUT 0xd0580000

.section .text
.global _start
.global bit_reverse_indices
.global bit_reverse_index
.global normalize_values
.global compute_twiddle_vector
.global FFT
.global printToLogVectorized
.global _finish

_start:
    la s0, cos_table  # cos table
    la s1, sin_table  # sin table
    add t1, zero, 0x500
    la a0, signal
    lw a1, size
    call printToLogVectorized
    la a1, signal
    lw a0, size
    call FFT
    mv a0, a1
    lw a1, size
    call printToLogVectorized
    j _finish

FFT:
    li t3, 2  # constant pair-wise size
    addi sp, sp, -12   # reserving stack space
    li a3, 0x600       # loading starting point for imaginary array so no memory conflicts
    sw a1,  0(sp)
    sw a3,  4(sp)
    sw ra,  8(sp)
    call bit_reverse_indices  # reversing array values based on its indices

perform_butterflies:
    bgt t3, a0, end_butterfly  # i > size, where i = 2,4,8,16....
    li s3, 0                   # j counter
    srli t5, t3, 1             # j/2  since array[j: j+i] will be divided in two parts: u's and v's
    slli a7, t5, 2             # (j/2) * 4 memory steps

perform_butteflies_over_sets:
    bge s3, a0, end_set_computation # j > size, where j = 0, j+i, j+2i....
    slli s2, s3, 2    # j*4
    add s4, t5, s3   # n/2+j -> address for second half(real)
    slli s4, t5, 2    # (n/2+j)*4
    add a1, a1 ,s2  # real[j] -> address for first half(real)
    add a3, a3, s2  # imag[j] -> address for first half(imaginary)
    add a2, a1, s4  # real[len/2 + j] -> address for second half(real)
    add a4, a3, s4  # imag[len/2 + j] -> address for second half(imaginary)
    mv t6, t5       # t5 -> i/2
    li s5, 0

compute_pair_over:
    vsetvli t0, t6, e32, m1
    sw a7, 12(sp)
    call compute_twiddle_vector # returns t1 = [cos(2*pi*k/i)] t2 = [sin(2*pi*k/i)] for i= 0,1,2,3...vlen 
    mv t0, t0
    lw a7, 0(a1)
    lw a7, 12(sp)
    vle32.v v1, (a1)  # real first half,  u_real
    vle32.v v2, (a2)  # real second half, v_real
    vle32.v v3, (a3)  # imag first half,  u_imag
    vle32.v v4, (a4)  # imag second half, v_imag
    vle32.v v5, (t1)  # twiddle factors(real part) vector    w_real
    vle32.v v6, (t2)  # twiddle factors(imag part) vector    w_imag
    vfmul.vv v7, v5, v2   # w_real * v_real
    vfmul.vv v8, v6, v4   # w_imag * v_imag 
    vfsub.vv v9, v7, v8   # w_real * v_real - w_imag * v_imag -> v*w real part
    vfmul.vv v7, v5, v4   # w_real * v_imag
    vfmul.vv v8, v6, v2   # w_imag * v_real 
    vfadd.vv v10, v7, v8   # w_real * v_imag - w_imag * v_real -> v*w imag part
    vfadd.vv v11, v1, v9  # (real)u + w*v
    vfadd.vv v12, v3, v10  #  (imag)u + w*v 
    vfsub.vv v13, v1, v9  # (real)u - w*v
    vfsub.vv v14, v3, v10  #  (imag)u - w*v
    vse32.v v11, (a1)    # saving first-half(real) butterfly 
    vse32.v v13, (a2)    # saving second-half(real) butterfly 
    vse32.v v12, (a3)    # saving first-half(imag) butterfly 
    vse32.v v14, (a4)    # saving second-half(real) butterfly 
    slli s4, t0, 2  # vlen * 4
  
 increment:
    add a1, a1, s4
    add a2, a2, s4
    add a3, a3, s4
    add a4, a4, s4
    sub t6, t6, t0    # counter to check how many sets elements left
    add s5, s5, t0    # k-counter to tell where to start twiddle calculation from
    bgt t6, zero, compute_pair_over  # checks if more of set is left
   
next_pair:
    lw a1,  0(sp)   # loading real array starting point
    lw a3,  4(sp)   # loading imaginary array starting point
    add s3, s3, t3  # j = j + i
    addi s5, zero, 0 # resets the k counter required for twiddle
    j perform_butteflies_over_sets   # perform butterfly on next set
   
end_set_computation:   # marks the end of i-th butterfly  
    slli t3, t3, 1       # i = i*2
    j perform_butterflies  # perform next butterfly

end_butterfly:        # marks the end of all butterflies
    call normalize_values  # normalize the real and imag values
    lw a1, 0(sp)     # restore pointers
    lw a3, 4(sp)
    lw ra, 8(sp)
    addi sp, sp, 12   # releasing stack space
    ret

# === Bit Reverse Indices ===
# reverses a array values according to bit-reverse-permutation
# Input: a1 has the array and a0 has the size
# Output: a1 has now all the values adjusted in array according to reverse indices

bit_reverse_indices:
    addi sp, sp, -52
    sw t4, 48(sp)
    sw t3, 44(sp)
    sw t2, 40(sp)
    sw t1, 36(sp)
    sw t0, 32(sp)
    sw a1, 28(sp)
    sw a0, 24(sp)
    sw ra, 20(sp)
    sw s0, 16(sp)
    sw s1, 12(sp)
    sw s2, 8(sp)
    sw s3, 4(sp)
    sw s4, 0(sp)

    mv s0, a1             # s0 = base address
    mv s1, a0             # s1 = n
    # Compute log2(n)
    mv t0, s1
    li s2, 0
    
log2_loop:
    srli t0, t0, 1
    beq t0, zero, log2_done
    addi s2, s2, 1
    j log2_loop
log2_done:

    li s3, 1              # i = 1

reverse_loop:
    bge s3, s1, reverse_done

    mv a0, s3
    mv a1, s2
    jal bit_reverse_index
    mv s4, a0             # j

    bge s3, s4, next_index

    # Swap array[i] and array[j]
    slli t0, s3, 2
    add t0, s0, t0
    lw t1, 0(t0)

    slli t2, s4, 2
    add t2, s0, t2
    lw t3, 0(t2)

    sw t3, 0(t0)
    sw t1, 0(t2)

next_index:
    addi s3, s3, 1
    j reverse_loop

reverse_done:
    lw s4, 0(sp)
    lw s3, 4(sp)
    lw s2, 8(sp)
    lw s1, 12(sp)
    lw s0, 16(sp)
    lw ra, 20(sp)
    lw a0, 24(sp)
    lw a1, 28(sp)
    lw t0, 32(sp)
    lw t1, 36(sp)
    lw t2, 40(sp)
    lw t3, 44(sp)
    lw t4, 48(sp)
    addi sp, sp, 52
    ret

# === Bit Reverse Helper ===
# Reverses a single numbers binary digits
# Input: a single number in t1 register
# Output: reversed binary version of that number in a0
bit_reverse_index:
    li t0, 0
    mv t1, a0
    li t2, 0

reverse_bit_loop:
    beq t2, a1, reverse_bit_done
    andi t4, t1, 1
    slli t0, t0, 1
    or t0, t0, t4
    srli t1, t1, 1
    addi t2, t2, 1
    j reverse_bit_loop

reverse_bit_done:
    mv a0, t0
    ret

# Function: normalize_values
# computes the normalized version of real and imaginary values as the formula output[[i] = sqrt( real[i]^2 + imag[i]^2 )
# Input: arrays for real and imaginary values are provided in a1 and a3
# Output: a normalized arrays of values given in a1

normalize_values:
     vsetvli t0, a0, e32, m1
     li t4, 0
     slli t5, t0, 2  # VLEN * 4  mem offset
normalize_loop:
     bge t4, a0, end_normalize
     
     vle32.v v1, (a1)    # load real values
     vle32.v v2, (a3)    # load imag values
     
     vfmul.vv v3, v1, v1  # (real[i])^2
     vfmul.vv v4, v2, v2  # (imag[i])^2
     
     vfadd.vv v5, v3, v4  # real^2 + imag^2
     
     vfsqrt.v v6, v5  # sqrt( real^2 + imag^2 )
      
     vse32.v v6, (a1)  # store the normalized vals
     
     add a1, a1, t5    # update pointers
     add a3, a3, t5    # update pointers
     
     add t4, t4, t0  # t4 = t4 + vlen  
     
     j normalize_loop
end_normalize:
   ret

# Function: compute_twiddle_vector
# computes the twiddle factors for FFT
compute_twiddle_vector:
     li s7, 1024
     li a7, 1023
     li s8, 4
     fcvt.s.w f1, t3  # i in float
     fcvt.s.w f2, s7  # 1024 in float 
     fcvt.s.w f3, s8
     vid.v     v2                 # v2[k]=k for k=0...vlen-1
     vadd.vx v2, v2, s5           # v2[k] = (k+s)
     vfcvt.f.xu.v v2, v2
     vfdiv.vf v2, v2, f1           # v2[k] = (k+s)/i for k = 0,1...vlen-1
     vfmul.vf v2, v2, f2           # v2[k] = (k+s/i)*1024 for k = 0,1...vlen-1
     vfcvt.xu.f.v v3, v2           # convert from floats to int
     mv s7, s9
     vmv.v.x   v4, a7               # Fill v1 with the value 1023
     vand.vv   v5, v4, v3            # v2[i] = v3[i] & 1023
     vfmul.vf v5, v5, f3          # multiply by 4 for byte segmentation
     
     vluxei32.v v6, (s0), v5           # load multiple cosine values
     vluxei32.v v7, (s1), v5           # load multiple sin values
     
     vse32.v v6, (t1)             # store at t0(real part)
     vse32.v v7, (t2)             # store at t1(imag part)
     ret

# Function: printToLogVectorized
# helps python script in reading data from the log file and performing operations on it
# Input: memory located at a0
# Output: None

printToLogVectorized:        
    addi sp, sp, -4
    sw a0, 0(sp)
    
    li t0, 0x123                 # Pattern for help in python script
    li t0, 0x456                 # Pattern for help in python script
    mv a1, a1                   # moving size to get it from log 
    li t0, 0		                # load i = 0
printloop:
    vsetvli t3, a1, e32           # Set VLEN based on a1
    slli t4, t3, 2                # Compute VLEN * 4 for address increment

    vle32.v v1, (a0)              # Load real[i] into v1
    add a0, a0, t4                # Increment pointer for real[] by VLEN * 4
    add t0, t0, t3                # Increment index

    bge t0, a1, endPrintLoop      # Exit loop if i >= size
    j printloop                   # Jump to start of loop
endPrintLoop:
    li t0, 0x123                    # Pattern for help in python script
    li t0, 0x456                    # Pattern for help in python script
    lw a0, 0(sp)
    addi sp, sp, 4
    li t0, 0xfeffeddc              # Pattern to help python script
    jr ra

# Function: _finish
# VeeR Related function which writes to to_host which stops the simulator
_finish:
    li x3, 0xd0580000
    addi x5, x0, 0xff
    sb x5, 0(x3)
    beq x0, x0, _finish

    .rept 100
        nop
    .endr

.data
## ALL DATA IS DEFINED HERE LIKE MATRIX, CONSTANTS ETC
## PLEASE DONT CHANGE ABOVE CODE ALL CHANGES TO BE DONR BELOW THIS POINT
## DATA DEFINE START

.equ mat_size, 1024

size: .word mat_size

# Sample input signal data: 128 floats
signal:
	.float  0.0000,  0.8778,  1.5174,  1.7596,  1.5759,  1.0729,  0.4497, -0.0760, -0.3464, -0.3182 
	.float -0.0732,  0.2186,  0.3678,  0.2433, -0.1733, -0.7748, -1.3635, -1.7209, -1.6883, -1.2270 
	.float -0.4370,  0.4725,  1.2528,  1.6981,  1.7134,  1.3430,  0.7495,  0.1520, -0.2542, -0.3665 
	.float -0.2085,  0.0852,  0.3244,  0.3414,  0.0591