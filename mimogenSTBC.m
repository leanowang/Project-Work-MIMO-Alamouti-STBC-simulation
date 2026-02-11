% Antonio Santos da Silva - 2025
% QPSK signal generation with preamble and export to VSG.

clear; clc;
message             = 'leano;';
sampleRate          = 25e6;          
symbolRate          = 2.5e6;           
samplesPerSymbol    = sampleRate / symbolRate;
rolloff             = 0.35;
filterSpan          = 6;               

pn = comm.PNSequence( ...
    'Polynomial',[8 6 5 4 1 0], ...
    'SamplesPerFrame',255, ...
    'InitialConditions',[1 0 0 0 0 0 0 0]);

preambleBits = pn();          % 255×1 列向量
preambleBits = preambleBits(:);   % 明确保证是列向量（好习惯）

if mod(length(preambleBits), 2) ~= 0
    preambleBits(end+1) = 0;  % 纵向补 1 个 bit
end

preambleBits = logical(preambleBits);

% messageBits = de2bi(uint8(message), 8, 'left-msb')';
% bits = messageBits(:);
% if mod(length(bits), 2) ~= 0
%     bits = [bits; 0];
% end
% bits = logical(bits);
% qpskMod = comm.QPSKModulator('BitInput', true, 'PhaseOffset', pi/4);
% preambleSymbols = qpskMod(preambleBits);
% messageSymbols = qpskMod(bits);
messageBits = de2bi(uint8(message), 8, 'left-msb')';
msgBits = logical(messageBits(:));

% ===== 长度头：消息字节数（uint16） =====
msgLenBytes = uint16(length(message));
lenBits = de2bi(msgLenBytes, 16, 'left-msb').';
lenBits = logical(lenBits(:));

% ===== payload bits = LEN + MSG =====
payloadBits = [lenBits; msgBits];
if mod(length(payloadBits),2) ~= 0
    payloadBits = [payloadBits; 0];
end

qpskMod = comm.QPSKModulator('BitInput', true, 'PhaseOffset', pi/4);
preambleSymbols = qpskMod(preambleBits);
messageSymbols  = qpskMod(payloadBits);   % 用 payloadBits 调制

Lp = length(preambleSymbols);
tx1_pre = [preambleSymbols; zeros(Lp,1)];  % 时隙A: Tx1发, Tx2静默
tx2_pre = [zeros(Lp,1);  preambleSymbols]; % 时隙B: Tx2发, Tx1静默

%  payload：only messageSymbols  Alamouti 
msg = messageSymbols;
if mod(length(msg),2)~=0
    msg = [msg; 0];
end
numPairs = length(msg)/2;
s1 = msg(1:2:end);
s2 = msg(2:2:end);

tx1_payload = zeros(2*numPairs,1);
tx2_payload = zeros(2*numPairs,1);
for k = 1:numPairs
    % timeslot 1
    tx1_payload(2*k-1) =  s1(k);
    tx2_payload(2*k-1) =  s2(k);
    % timeslot 2
    tx1_payload(2*k)   = -conj(s2(k));
    tx2_payload(2*k)   =  conj(s1(k));
end

% combine pre and payload
tx1_syms = [tx1_pre; tx1_payload];
tx2_syms = [tx2_pre; tx2_payload];

% ===== 在进入 Tx RRC 前补零符号：flush + 留足够尾巴 =====
padSyms = 2*filterSpan;          % 建议至少 12 个符号（很安全）
tx1_syms_f = [tx1_syms; zeros(padSyms,1)];
tx2_syms_f = [tx2_syms; zeros(padSyms,1)];

rrcFilter = comm.RaisedCosineTransmitFilter( ...
    'RolloffFactor', rolloff, ...
    'FilterSpanInSymbols', filterSpan, ...
    'OutputSamplesPerSymbol', samplesPerSymbol);

tx1Signal = rrcFilter(tx1_syms_f);
release(rrcFilter);
tx2Signal = rrcFilter(tx2_syms_f);

tx1Signal = tx1Signal./max(abs(tx1Signal));
tx2Signal = tx2Signal./max(abs(tx2Signal));

% 再归一化（可选，放前放后都行，只要两路一致）
tx1Signal = tx1Signal ./ max(abs(tx1Signal));
tx2Signal = tx2Signal ./ max(abs(tx2Signal));
%CSV
iData1 = real(tx1Signal); qData1 = imag(tx1Signal);
iData2 = real(tx2Signal); qData2 = imag(tx2Signal);
writematrix([iData1 qData1], 'leano1.csv');
writematrix([iData2 qData2], 'leano2.csv');

disp('已保存 iq_Tx1.csv 和 iq_Tx2.csv（preamble未编码，payload为Alamouti）');