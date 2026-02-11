%% rx_full_decode_with_cfo.m
% 一键：稳健读取 IQ -> DC+归一化 -> RRC -> PN检测 -> CFO估计补偿 ->
% 两段preamble估计h1/h2 -> Alamouti合并 -> QPSK 4旋转尝试 -> LEN+MSG解析

clear; clc; close all;

%% ============ 1) 你的接收文件 ============
rxFile = '20260209.txt';   % <-- 改成你的文件名（.txt / .csv 都行）

%% ============ 2) 系统参数（按你的TX/采样设置改） ============
sampleRate       = 25e6;
symbolRate       = 2.5e6;
samplesPerSymbol = sampleRate / symbolRate;  % 10（确保是整数）
if abs(samplesPerSymbol - round(samplesPerSymbol)) > 1e-9
    error('samplesPerSymbol is not integer. Check sampleRate/symbolRate.');
end
samplesPerSymbol = round(samplesPerSymbol);

phaseOffset      = pi/4;   % 你的QPSK调制 PhaseOffset
rolloff          = 0.35;
filterSpan       = 6;      % RRC span in symbols

doRRC = true;

%% ============ 3) 读 IQ（最关键） ============
[I, Q] = readIQ_robust(rxFile);
rxSignal = complex(I, Q);

fprintf('Loaded file: %s\n', rxFile);
fprintf('N=%d samples\n', numel(rxSignal));
fprintf('std(I)=%.6g, std(Q)=%.6g, ratio std(I)/std(Q)=%.6g\n', std(I), std(Q), std(I)/(std(Q)+eps));
fprintf('mean(I)=%.6g, mean(Q)=%.6g\n', mean(I), mean(Q));
fprintf('min/max I: [%.6g, %.6g]\n', min(I), max(I));
fprintf('min/max Q: [%.6g, %.6g]\n', min(Q), max(Q));
fprintf('first 10 rows (I Q):\n');
disp([I(1:min(10,end)) Q(1:min(10,end))]);

%% ============ 4) 原始IQ散点 ============
figure('Name','IQ scatter - raw');
plot(real(rxSignal), imag(rxSignal), '.', 'MarkerSize', 4);
grid on; axis equal;
xlabel('I'); ylabel('Q');
title('IQ scatter (RAW)'); 

%% ============ 5) 去 DC + 归一化 ============
rx0 = rxSignal - mean(rxSignal);
rx0 = rx0 ./ (max(abs(rx0)) + 1e-12);

figure('Name','IQ scatter - DC removed + normalized');
plot(real(rx0), imag(rx0), '.', 'MarkerSize', 4);
grid on; axis equal;
xlabel('I'); ylabel('Q');
title('IQ scatter (DC removed + normalized)');

%% ============ 6) （可选）功率随时间 ============
p = abs(rx0).^2;
p_s = movmean(p, 2000);
figure('Name','Smoothed power');
plot(p_s); grid on;
xlabel('sample index'); ylabel('|x|^2 (smoothed)');
title('Smoothed power');

%% ============ 7) RRC 滤波 ============
if doRRC
    rrcRx = comm.RaisedCosineReceiveFilter( ...
        'RolloffFactor', rolloff, ...
        'FilterSpanInSymbols', filterSpan, ...
        'InputSamplesPerSymbol', samplesPerSymbol, ...
        'DecimationFactor', 1);

    y = rrcRx(rx0);
else
    y = rx0;
end

% 群时延（样点）
gd = filterSpan * samplesPerSymbol;

%% ============ 8) 生成 PN preamble（务必与TX一致） ============
% ---- Tx side PN ----
pn = comm.PNSequence( ...
    'Polynomial',[8 6 5 4 1 0], ...
    'SamplesPerFrame',255, ...
    'InitialConditions',[1 0 0 0 0 0 0 0]);

preambleBits = pn();
preambleBits = preambleBits(:);
if mod(length(preambleBits),2)~=0
    preambleBits(end+1)=0;
end
preambleBits = logical(preambleBits);

qpskMod = comm.QPSKModulator('BitInput',true,'PhaseOffset',phaseOffset);
preambleSyms = qpskMod(preambleBits);
Lp = length(preambleSyms);

fprintf('PN preamble: bits=%d, Lp=%d symbols\n', length(preambleBits), Lp);

%% ============ 9) PN preamble detection（扫描k相位） ============
best = struct('k',0,'n0',1,'peak',-Inf,'match',-Inf,'phaseDeg',0);
bestCorrAbs = [];
bestCorrC = [];

for k = 0:samplesPerSymbol-1
    rxSyms_k = y(gd+1+k : samplesPerSymbol : end);
    if numel(rxSyms_k) < Lp + 10
        continue;
    end

    tpl = flipud(conj(preambleSyms(:)));
    num = conv(rxSyms_k, tpl, 'valid');
    rxPow = conv(abs(rxSyms_k).^2, ones(Lp,1), 'valid');
    tplPow = sum(abs(preambleSyms).^2) + 1e-12;

    corrC   = num ./ sqrt(rxPow*tplPow + 1e-12);
    corrAbs = abs(corrC);

    [pk, idx] = max(corrAbs);

    % pnMatch：峰值处去复增益后，与参考相似度(0..1)
    seg = rxSyms_k(idx:idx+Lp-1);
    g = (seg' * preambleSyms) / (preambleSyms' * preambleSyms + 1e-12);
    seg2 = seg / (g + 1e-12);
    pnMatch = abs(seg2' * preambleSyms) / (norm(seg2)*norm(preambleSyms) + 1e-12);

    if pk > best.peak
        best.k = k;
        best.n0 = idx;
        best.peak = pk;
        best.match = pnMatch;
        best.phaseDeg = rad2deg(angle(corrC(idx)));
        bestCorrAbs = corrAbs;
        bestCorrC = corrC;
    end
end

fprintf('PN search result: best_k=%d, bestStart=%d (symbol index), |corr|=%.4f, pnMatch=%.4f, phase=%.1f deg\n', ...
    best.k, best.n0, best.peak, best.match, best.phaseDeg);

figure('Name','PN normalized correlation |corr|');
plot(bestCorrAbs); grid on;
xlabel('symbol index n'); ylabel('|corr(n)|');
title(sprintf('PN corr |corr| (best k=%d, peak=%.3f, pnMatch=%.3f)', best.k, best.peak, best.match));
hold on; plot(best.n0, bestCorrAbs(best.n0), 'ro'); hold off;

% 叠加图：检测到的preamble段 vs 参考
rxSyms_best = y(gd+1+best.k : samplesPerSymbol : end);
seg = rxSyms_best(best.n0:best.n0+Lp-1);
g = (seg' * preambleSyms) / (preambleSyms' * preambleSyms + 1e-12);
seg2 = seg / (g + 1e-12);

figure('Name','Preamble overlay (detected segment vs PN)');
plot(real(seg2), imag(seg2), '.'); hold on;
plot(real(preambleSyms), imag(preambleSyms), '.');
grid on; axis equal;
legend('rx segment (gain removed)','PN reference');
title('Detected preamble segment vs PN reference');

%% ============ 10) CFO 估计 + 补偿（关键新增） ============
% 用 e[n] = (seg去增益) * conj(preamble) 拿到相位斜坡，线性拟合斜率 => CFO
e = seg2 .* conj(preambleSyms);
ph = unwrap(angle(e));
n  = (0:Lp-1).';
pp = polyfit(n, ph, 1);
omega_sym = pp(1);                       % rad / symbol
cfo_hz = omega_sym/(2*pi) * symbolRate;  % Hz

fprintf('Estimated CFO ≈ %.3f Hz (symbol-domain, from PN)\n', cfo_hz);

% 对整个符号流做补偿（在RRC+抽样之后的符号域）
Nsym = length(rxSyms_best);
rot = exp(-1j * omega_sym * (0:Nsym-1).');
rxSyms_best_cfo = rxSyms_best .* rot;

% 可视化：补偿后preamble是否更"站住"
seg_cfo = rxSyms_best_cfo(best.n0:best.n0+Lp-1);
g2 = (seg_cfo' * preambleSyms) / (preambleSyms' * preambleSyms + 1e-12);
seg2_cfo = seg_cfo / (g2 + 1e-12);

figure('Name','Preamble overlay AFTER CFO comp');
plot(real(seg2_cfo), imag(seg2_cfo), '.'); hold on;
plot(real(preambleSyms), imag(preambleSyms), '.');
grid on; axis equal;
legend('rx segment (gain removed, CFO comp)','PN reference');
title('Preamble overlay after CFO compensation');

%% ============ 11) 用 best_k/bestStart 解一帧：两段preamble估计h1/h2 ============
best_k     = best.k;
frameStart = best.n0;

% 注意：后续全部用 CFO 补偿后的 rxSyms
rxSyms = rxSyms_best_cfo;

% 两个preamble连续：A then B
idxA = frameStart;
idxB = frameStart + Lp;

if idxB + Lp - 1 > length(rxSyms)
    error('Not enough symbols for two preambles at detected start.');
end

yA = rxSyms(idxA : idxA+Lp-1);
yB = rxSyms(idxB : idxB+Lp-1);

den = sum(abs(preambleSyms).^2) + 1e-12;
h1 = sum(yA .* conj(preambleSyms)) / den;
h2 = sum(yB .* conj(preambleSyms)) / den;

fprintf('h1=%.4g%+.4gi  |h1|^2=%.4g\n', real(h1), imag(h1), abs(h1)^2);
fprintf('h2=%.4g%+.4gi  |h2|^2=%.4g\n', real(h2), imag(h2), abs(h2)^2);

%% ============ 12) payload 截取 + Alamouti 合并 ============
payloadStart = frameStart + 2*Lp;

winSyms = 60;   % 先取长一点；最好你后续按帧长精确切
if payloadStart > length(rxSyms)
    error('payloadStart beyond rxSyms length.');
end
payloadSyms = rxSyms(payloadStart : min(end, payloadStart+winSyms-1));

% Alamouti needs even length
if mod(length(payloadSyms),2)~=0
    payloadSyms = payloadSyms(1:end-1);
end

numBlocks = length(payloadSyms)/2;
Hn = abs(h1)^2 + abs(h2)^2 + 1e-12;

s1_hat = zeros(numBlocks,1);
s2_hat = zeros(numBlocks,1);

for kk = 1:numBlocks
    y1 = payloadSyms(2*kk-1);
    y2 = payloadSyms(2*kk);

    s1_hat(kk) = (conj(h1)*y1 + h2*conj(y2)) / Hn;
    s2_hat(kk) = (conj(h2)*y1 - h1*conj(y2)) / Hn;
end

rxHat = reshape([s1_hat.'; s2_hat.'], [], 1);

figure('Name','Payload constellation after Alamouti (CFO comp)');
plot(real(rxHat(1:min(5000,end))), imag(rxHat(1:min(5000,end))), '.', 'MarkerSize', 4);
grid on; axis equal;
xlabel('I'); ylabel('Q');
title('Payload constellation after Alamouti combining (after CFO comp)');

%% ============ 13) QPSK demod with 4 rotations, then parse LEN+MSG ============
qpskDemod = comm.QPSKDemodulator('BitOutput', true, 'PhaseOffset', phaseOffset);

bestMsg = "";
bestRot = -1;
bestLen = NaN;

for rotk = 0:3
    rxRot = rxHat * exp(-1j*rotk*pi/2);
    bits = qpskDemod(rxRot);

    if length(bits) < 16
        continue;
    end

    lenBits = bits(1:16).';
    msgLenBytes = bi2de(lenBits,'left-msb');
    needBits = 16 + double(msgLenBytes)*8;

    % 合理性过滤
    if msgLenBytes==0 || msgLenBytes > 200
        continue;
    end
    if length(bits) < needBits
        continue;
    end

    msgBits = bits(17 : 16+double(msgLenBytes)*8);
    msgBits = msgBits(1:floor(length(msgBits)/8)*8);
    B = reshape(msgBits, 8, []).';
    bytes = bi2de(B,'left-msb');
    msg = char(bytes).';

    bestMsg = msg;
    bestRot = rotk;
    bestLen = msgLenBytes;
    break;
end

fprintf('Decoded: rot=%d, LEN=%d, msg="%s"\n', bestRot, bestLen, bestMsg);

%% ======================= 本文件末尾：函数 =======================
function [I, Q] = readIQ_robust(fname)
% Robustly read 2-column IQ from txt/csv with:
% - comma/semicolon/space delimiters
% - decimal comma (e.g., 0,00123) formats
% - optional header / extra text lines

    raw = fileread(fname);
    raw = strrep(raw, sprintf('\r'), ''); % normalize CRLF

    % Detect decimal comma: digit,comma,digit
    isDecComma = ~isempty(regexp(raw, '\d,\d', 'once'));

    if isDecComma
        % Convert decimal comma to decimal dot: -12,345 -> -12.345
        raw = regexprep(raw, '(-?\d+),(\d+)', '$1.$2');
    end

    % Unify delimiters to spaces
    raw = strrep(raw, ',', ' ');
    raw = strrep(raw, ';', ' ');
    raw = strrep(raw, sprintf('\t'), ' ');

    % Keep only numeric-relevant chars
    raw = regexprep(raw, '[^0-9eE\+\-\. \n]', ' ');

    % Parse all floats
    nums = sscanf(raw, '%f');

    if isempty(nums)
        error('No numeric data parsed. Check file encoding/content.');
    end
    if mod(numel(nums), 2) ~= 0
        error('Parsed an odd number of floats (%d). Not a 2-column IQ file.', numel(nums));
    end

    nums = reshape(nums, 2, []).';
    I = nums(:,1);
    Q = nums(:,2);

    if std(I)==0 && std(Q)>0
        warning('I is constant but Q varies. Parsing may still be wrong (column swap or format issue).');
    end
end
