/*
 * Copyright (c) 2024 Your Name
 * SPDX-License-Identifier: Apache-2.0
 */

`default_nettype none
module ERYTHwrapper (
    input  wire [7:0] ui_in,    // Dedicated inputs//I1
    output wire [7:0] uo_out,   // Dedicated outputs//OUTPUT
    input  wire [7:0] uio_in,   // IOs: Input path//I2
    output wire [7:0] uio_out,  // IOs: Output path
    output wire [7:0] uio_oe,   // IOs: Enable path (active high: 0=input, 1=output)
    input  wire       ena,      // always 1 when the design is powered, so you can ignore it
    input  wire       clk,      // clock
    input  wire       rst_n     // reset_n - low to reset
);

  // All output pins must be assigned. If not used, assign to 0.
  //assign uo_out  = ui_in + uio_in;  // Example: ou_out is the sum of ui_in and uio_in
  assign uio_out = 0;
  assign uio_oe  = 0;

  // List all unused inputs to prevent warnings
//  wire _unused = &{ena, clk, rst_n, 1'b0};

//endmodule
//module erythcrypt_final(
 //   input [7:0] I1,
  //  input [7:0] I2,
//    input CLK,
 //   input Reset,
 //   input [3:0] Control,
 //   output reg [7:0] OUTPUT
//);

    wire [7:0] Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10, Y11, Y12;
    wire [7:0] I3, I4, I5;
   
    assign I3 = 8'b00000011;
    assign I4 = 8'b00000101;
    assign I5 = 8'b00011001;
   
    gcd2 g1(.a(ui_in), .b(uio_in), .clk(clk), .reset(rst_n), .gcd2(Y1));
    mod3 m1(.a(ui_in), .b(uio_in), .clk(clk), .rst(rst_n), .rem(Y2));
    rightshift l1(.A(ui_in), .S(uio_in[0]), .Y(Y3));
    rightshift1 r1(.A(ui_in), .S(uio_in[0]), .Y(Y4));
    matrixadd ma1(.A(ui_in), .B(uio_in), .Res(Y5));
    matrixmult mm1(.A(ui_in), .B(uio_in), .Res(Y6));
    mul_inv mi1(.a(ui_in), .b(uio_in), .n(I3), .clk(clk), .reset(rst_n), .mul_in(Y7));
    congruence c1(.a(ui_in), .b(uio_in), .n(I4), .clk(clk), .reset(rst_n), .congruent(Y8));
    mod_inverse md1(.a(ui_in), .m(uio_in), .inverse(Y9));
    exp e1(.a(ui_in), .b(uio_in), .m(I5), .clk(clk), .out(Y10));
    flt f1(.a(ui_in), .p(uio_in), .clk(clk), .reset(rst_n), .result(Y11));
    //.done(Y11));
//     Instantiate the 16-to-1 multiplexer
        wire [7:0] mux_output;
        mux_16to1 mux (
            .A(Y1), .B(Y2), .C(Y3), .D(Y4), .E(Y5), .F(Y6), .G(Y7), .H(Y8),
            .I(Y9), .J(Y10), .K(Y11),
//             .M(Y13), .N(Y14), .O(Y15), .P(Y16),
            .So(Control),
            .Y(mux_output)
        );
       
        always @(posedge clk) begin
            if (rst_n==0)
                OUTPUT <= 0;
            else
                OUTPUT <= mux_output;
        end
    endmodule
   




module modulo(
    input [7:0] a,     // Dividend (numerator)
    input [7:0] b,     // Divisor (denominator)
    output reg [7:0] result // Remainder
);

    always @(*) begin
        if (b != 0) begin
            result = a % b;  // Calculate a modulo b
        end else begin
            result = 0;      // If divisor is 0, result is set to 0 (undefined behavior for division by 0)
        end
    end
endmodule

module comparator(
    input [7:0] a,     // 8-bit input a
    input [7:0] b,     // 8-bit input b
    output reg a_lt_b, // Output: a < b
    output reg a_eq_b, // Output: a == b
    output reg a_gt_b  // Output: a > b
);

    always @(*) begin
        if (a < b) begin
            a_lt_b = 1;
            a_eq_b = 0;
            a_gt_b = 0;
        end else if (a == b) begin
            a_lt_b = 0;
            a_eq_b = 1;
            a_gt_b = 0;
        end else begin
            a_lt_b = 0;
            a_eq_b = 0;
            a_gt_b = 1;
        end
    end

endmodule


module gcd2 (
    input [7:0] a, b,        // 8-bit inputs a and b
    input clk,               // Clock input
    input reset,             // Reset signal
    output reg [7:0] gcd2    // 8-bit output for GCD
);

    reg [7:0] min;           // Register to store minimum of a and b
    reg [7:0] counter;       // Counter register
    wire a_lt_b, a_eq_b, a_gt_b;  // Comparator outputs
    wire [7:0] mod_a, mod_b;      // Modulo results

    // Instantiate the comparator to compare a and b
    comparator comp (
        .a(a),
        .b(b),
        .a_lt_b(a_lt_b),
        .a_eq_b(a_eq_b),
        .a_gt_b(a_gt_b)
    );

    // Instantiate mod_op_sub for modulo operations a % counter and b % counter
    modulo moda (
        .a(a),
        .b(counter),
        .result(mod_a)
    );

    modulo modb (
        .a(b),
        .b(counter),
        .result(mod_b)
    );

    always @(posedge clk) begin
        if (reset==0) begin
            gcd2 <= 8'd0;       // Reset GCD output to 0
            min <= 8'd0;        // Reset min
            counter <= 8'd1;    // Initialize counter to 1
        end else begin
            // Select min based on comparison using comparator module
            if (a_lt_b) begin
                min <= a;  // If a < b, min = a
            end else if (a_gt_b) begin
                min <= b;  // If a > b, min = b
            end

            // Counter to iterate from 1 to min
            if (counter <= min) begin
                // If both a % counter == 0 and b % counter == 0, gcd1 = counter
                if ((mod_a == 8'd0) && (mod_b == 8'd0)) begin
                    gcd2 <= counter;
                end
                counter <= counter + 1;  // Increment counter
            end
        end
    end
endmodule


/////////////////////////////////////////////////////////////////////////////////////////////////////MODULO
module mod3(
    input [7:0] a,
    input [7:0] b,
    input clk,
    input rst,
    output reg [7:0] rem
);

    always @(posedge clk) begin
        if (rst)
            rem <= 8'b0; // Ensure remainder is also reset
         else begin
            if(a>b)
               rem<=a%b;
            else
               rem<=b%a;
         end
      end
endmodule

////////////////////////////////////////////////////////////////////////////////////////////////////LEFTSHIFT
module rightshift(
input [7:0]A,input S,
output [7:0]Y
    );
    wire So;
    wire B;
    mux_1 m9(A[0],1'b0,S,Y[0]);
    mux_1 m10(A[1],A[0],S,Y[1]);
    mux_1 m11(A[2],A[1],S,Y[2]);
    mux_1 m12(A[3],A[2],S,Y[3]);
    mux_1 m13(A[4],A[3],S,Y[4]);
    mux_1 m14(A[5],A[4],S,Y[5]);
    mux_1 m15(A[6],A[5],S,Y[6]);
    mux_1 m16(A[7],A[6],S,Y[7]);
endmodule

/////////////////////////////////////////////////////////////////////////////////////////////////////RIGHTSHIFT
module rightshift1(
input [7:0]A,input S,
output [7:0]Y
    );
    wire So;
    wire B;
    mux_1 m1(A[0],A[1],S,Y[0]);
    mux_1 m2(A[1],A[2],S,Y[1]);
    mux_1 m3(A[2],A[3],S,Y[2]);
    mux_1 m4(A[3],A[4],S,Y[3]);
    mux_1 m5(A[4],A[5],S,Y[4]);
    mux_1 m6(A[5],A[6],S,Y[5]);
    mux_1 m7(A[6],A[7],S,Y[6]);
    mux_1 m8(A[7],1'b0,S,Y[7]);
endmodule

//////////////////////////////////////////////////////////////////////////////////////////////////////MUX_1
module mux_1(
    input A,B,So,
    output Y
    );
    assign Y=(~So&A)|(So&B);
endmodule

//////////////////////////////////////////////////////////////////////////////////////////////////////MATRIX ADD
module matrixadd(A,B,Res );
input [7:0] A;
input [7:0] B;
output [7:0] Res;

reg [7:0] Res;
reg [7:0] A1 [0:1][0:1];
reg [7:0] B1 [0:1][0:1];
reg [7:0] Res1 [0:1][0:1];
integer i,j;

always@ (A or B)
begin

{A1[0][0],A1[0][1],A1[1][0],A1[1][1]} = A;
{B1[0][0],B1[0][1],B1[1][0],B1[1][1]} = B;
i = 0;
j = 0;

{Res1[0][0],Res1[0][1],Res1[1][0],Res1[1][1]} = 8'd0;

for(i=0;i<2;i=i+1)
for(j=0;j<2;j=j+1)

Res1[i][j] = (A1[i][j] + B1[i][j]);

Res= {Res1[0][0],Res1[0][1],Res1[1][0],Res1[1][1]};

end
endmodule

/////////////////////////////////////////////////////////////////////////////////////////////////////////////MATRIX MULTIPLICATION
module matrixmult(A,B,Res );
input [7:0] A;
input [7:0] B;
output [7:0] Res;

//internal variables
reg [7:0] Res;
reg [7:0] A1 [0:1][0:1];


reg [7:0] B1 [0:1][0:1];
reg [7:0] Res1 [0:1][0:1];
integer i,j,k;

always@ (A or B)
begin
//Initialize the matrices-convert 1 D to 3D arrays
{A1[0][0],A1[0][1],A1[1][0],A1[1][1]} = A;
{B1[0][0],B1[0][1],B1[1][0],B1[1][1]} = B;
i = 0;
j = 0;
k = 0;
{Res1[0][0],Res1[0][1],Res1[1][0],Res1[1][1]} = 8'd0; //initialize to
//zeros.
//Matrix multiplication
for(i=0;i<2;i=i+1)
for(j=0;j<2;j=j+1)
for(k=0;k<2;k=k+1)
Res1[i][j] = Res1[i][j] + (A1[i][k] * B1[k][j]);
//final output assignment - 3D array to 1D array conversion.

Res = {Res1[0][0],Res1[0][1],Res1[1][0],Res1[1][1]};
end

endmodule

//////////////////////////////////////////////////////////////////////////////////////////////////MULTIPLICATIVE INVERSE
module mul_inv (
    input [7:0] a,
    input [7:0] b,
    input [7:0] n,
    input clk,
    input reset,
    output reg [7:0]mul_in
);
    wire [7:0] temp = 8'd1; // For 1 mod n operation
    wire [7:0] mult_out;
    wire [7:0]congruent;

    // Instantiate the multiplier module
    multiplier m1 (
        .a(a),
        .b(b),
        .multiply(mult_out),
        .clk(clk),
        .reset(reset)
    );

    // Instantiate the congruence module
    congruence c1 (
        .a(mult_out),
        .b(temp),
        .n(n),
        .clk(clk),
        .reset(reset),
        .congruent(congruent)
    );

    always @(posedge clk) begin
        if (reset==0) begin
            mul_in <= 0;
        end else begin
            mul_in <= congruent;
        end
    end
endmodule

//////////////////////////////////////////////////////////////////////////////////////////////////////////MULTIPLIER
module multiplier(a,b,multiply,clk,reset);
    input [7:0]a,b;
    input clk,reset;
    output reg [7:0] multiply;
    //initial begin multiply<=0; end
    always @(posedge clk)
        begin
            if(reset==0)
            begin
            multiply<=a*b;
            end
            else
            begin
             multiply<=0;
            end
        end
endmodule

////////////////////////////////////////////////////////////////////////////////////////////////////////////CONGRUENCE
module congruence (
    input [7:0] a,    // Input number a
    input [7:0] b,    // Input number b
    input [7:0] n,    // Modulus n
    input clk,        // Clock signal
    input reset,      // Reset signal
    output reg [7:0] congruent // Output signal indicating if a ? b (mod n)
);

    reg [7:0] a_mod_n, b_mod_n;

    always @(posedge clk) begin
        if (reset==0) begin
            congruent <= 0;
        end else begin
            // Compute the remainder of a and b when divided by n
            a_mod_n <= a % n;
            b_mod_n <= b % n;
           
            // Check if remainders are equal
            if (a_mod_n == b_mod_n) begin
                congruent <= 1;
            end else begin
                congruent <= 0;
            end
        end
    end

endmodule

/////////////////////////////////////////////////////////////////////////////////////////////////////MODULAR INVERSE
module mod_inverse(
    input [7:0] a,       // Input value
    input [7:0] m,       // Modulus value
    output reg[7:0] inverse // Output: Modular inverse of 'a' modulo 'm'
);

    reg [7:0] x, y;
    reg [7:0] x0, x1, y0, y1;
    reg [7:0] quotient, remainder, temp;
    reg [7:0] t;

    always @(*) begin
        x0 = 0; x1 = 1; y0 = 1; y1 = 0;
        x = a; y = m;
        if (y != 0) begin
            quotient = x / y;
            remainder = x % y;
            x = y;
            y = remainder;
            t = x0 - quotient * x1;
            x0 = x1;
            x1 = t;
            t = y0 - quotient * y1;
            y0 = y1;
            y1 = t;
        end
        inverse = x0 >= 0 ? x0 : x0 + m;
    end
endmodule

//////////////////////////////////////////////////////////////////////////////////////////////////MODULAR EXPONENTIATION
module exp(
    input [7:0]a,
    input [7:0]b,
    input [7:0]m,
    output reg [7:0]out,
    input clk
   
    );
reg [7:0]x1;
reg [7:0]y1;
always @(posedge clk) begin
       x1<=a-m;
       y1<=x1^b;
       out<=y1%m;
end
endmodule



module flt(
    input [7:0] a,
    input [7:0] p,
    input clk,
    input reset,
    output reg [7:0]result
    );
    reg[7:0]b;
    reg[7:0]c;
    reg[7:0]d;
    always@(posedge clk)begin
        if(reset==0)begin
           b<=0;
       end else begin
           if((a%p)!=0)begin
              b<=p-1;
              c<=a^b;
              d<=c%p;
           if(d==1)begin
              result<=1;
           end else begin
              result<=0;
           end
           end
          else
             result<=0;
        end
    end
endmodule



