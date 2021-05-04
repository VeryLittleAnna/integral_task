

section .data
    c11 dt 0.6
    c12 dd 3
    c21 dd 2
    c22 dd 1
    c31 dd -3
    fstr db `%Lf\n`, 0

section .text
global f1
global f1_derivative
global f2
global f2_derivative
global f3
global f3_derivative
extern printf


f1: ; f1 = 0.6x + 3
    push ebp
    mov ebp, esp
    finit
    fld tword[ebp + 8] ; x
    fld tword[c11]
    fmulp
    fild dword[c12]
    faddp
    mov esp, ebp
    pop ebp
    ret
    
f1_derivative: ; f1_d = 0.6
    push ebp
    mov ebp, esp
    finit
    fld tword[c11] 
    mov esp, ebp
    pop ebp
    ret
    
f2: ; f2 = (x - 2) ^ 3 - 1
    push ebp
    mov ebp, esp
    sub esp, 12
    finit 
    fld tword[ebp + 8]
    fild dword[c21] ;2
    fsubp
    fstp tword[ebp - 16] ; no fst for tword
    fld tword[ebp - 16]
    fld tword[ebp - 16]
    fld tword[ebp - 16]
    fmulp
    fmulp
    fild dword[c22] ; 1
    fsubp
    mov esp, ebp
    pop ebp
    ret

    


f2_derivative: ; f2_d = 3 * (x - 2) ^ 2
    push ebp
    mov ebp, esp
    sub esp, 12
    finit 
    fld tword[ebp + 8]
    fild dword[c21] ;2
    fsubp
    fstp tword[ebp - 16]
    fld tword[ebp - 16]
    fld tword[ebp - 16]
    fmulp
    fild dword[c12] ; 3
    fmulp 
    mov esp, ebp
    pop ebp
    ret
    
f3:
    push ebp
    mov ebp, esp
    finit 
    fild dword[c12] ; 3
    fld tword[ebp + 8]
    fdivp
    mov esp, ebp
    pop ebp
    ret
    
f3_derivative:
    push ebp
    mov ebp, esp
    finit 
    fild dword[c31] ; -3
    fld tword[ebp + 8]
    fdiv
    fdivp
    mov esp, ebp
    pop ebp
    ret