module div #(
    parameter WIDTH=32,  // width of numbers in bits (integer and fractional)
    parameter FBITS=29   // fractional bits within WIDTH
    ) (
    input wire clk,    // clock
    input wire rst,    // reset
    input wire start,  // start calculation
    output reg busy,   // calculation in progress
    output reg done,   // calculation is complete (high for one tick)
    output reg valid,  // result is valid
    output reg dbz,    // divide by zero
    output reg ovf,    // overflow
    input wire signed [WIDTH-1:0] a,   // dividend (numerator)
    input wire signed [WIDTH-1:0] b,   // divisor (denominator)
    output reg signed [WIDTH-1:0] val  // result value: quotient
    );

    localparam WIDTHU = WIDTH - 1;                 // unsigned widths are 1 bit narrower
    localparam FBITSW = (FBITS == 0) ? 1 : FBITS;  // avoid negative vector width when FBITS=0
    localparam SMALLEST = {1'b1, {WIDTHU{1'b0}}};  // smallest negative number
    localparam ITER = WIDTHU + FBITS;              // iteration count: unsigned input width + fractional bits

    localparam IDLE = 3'b000,
               INIT = 3'b001,
               CALC = 3'b010,
               ROUND = 3'b011,
               SIGN = 3'b100;

    reg [2:0] state;

    reg [$clog2(ITER):0] i;            // iteration counter (allow ITER+1 iterations for rounding)

    reg a_sig, b_sig, sig_diff;        // signs of inputs and whether different
    reg [WIDTHU-1:0] au, bu;           // absolute version of inputs (unsigned)
    reg [WIDTHU-1:0] quo, quo_next;    // intermediate quotients (unsigned)
    reg [WIDTHU:0] acc, acc_next;      // accumulator (unsigned but 1 bit wider)

    // input signs
    always @(*) begin
        a_sig = a[WIDTH-1];
        b_sig = b[WIDTH-1];
    end

    // Convert inputs to absolute values
    always @(*) begin
        au = (a_sig) ? -a : a;
        bu = (b_sig) ? -b : b;
    end

    // division algorithm iteration
    always @(*) begin
        if (acc >= {1'b0, bu}) begin
            acc_next = acc - bu;
            {acc_next[WIDTHU:0], quo_next} = {acc_next[WIDTHU-1:0], quo, 1'b1};
        end else begin
            {acc_next[WIDTHU:0], quo_next} = {acc, quo} << 1;
        end
    end

    // calculation state machine
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            state <= IDLE;
            busy <= 0;
            done <= 0;
            valid <= 0;
            dbz <= 0;
            ovf <= 0;
            val <= 0;
        end else begin
            done <= 0;
            case (state)
                IDLE: begin
                    if (start) begin
                        valid <= 0;
                        if (b == 0) begin  // divide by zero
                            state <= IDLE;
                            busy <= 0;
                            done <= 1;
                            dbz <= 1;
                            ovf <= 0;
                        end else begin
                            state <= INIT;
                            busy <= 1;
                            dbz <= 0;
                            ovf <= 0;
                        end
                    end
                end

                INIT: begin
                    state <= CALC;
                    ovf <= 0;
                    i <= 0;
                    {acc, quo} <= {{WIDTHU{1'b0}}, au};  // initialize calculation
                end

                CALC: begin
                    if (i == ITER-1) begin
                        if (acc_next[WIDTHU:0] > {1'b0, bu}) begin  // Check for overflow
                            ovf <= 1;
                            state <= IDLE;
                        end else begin
                            state <= ROUND;
                        end
                    end else begin
                        i <= i + 1;
                        acc <= acc_next;
                        quo <= quo_next;
                    end
                end

                ROUND: begin  // Gaussian rounding
                    state <= SIGN;
                    if (quo_next[0] == 1'b1) begin  // next digit is 1, so consider rounding
                        // round up if quotient is odd or remainder is non-zero
                        if (quo[0] == 1'b1 || acc_next[WIDTHU:1] != 0) quo <= quo + 1;
                    end
                end

                SIGN: begin  // adjust quotient sign if non-zero and input signs differ
                    state <= IDLE;
                    if (quo != 0) begin
                        val <= (a_sig ^ b_sig) ? -quo : quo;
                    end else begin
                        val <= 0;
                    end
                    busy <= 0;
                    done <= 1;
                    valid <= 1;
                end

                default: state <= IDLE;
            endcase
        end
    end
endmodule
