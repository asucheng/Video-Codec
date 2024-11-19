function reference_frames = A2_Q1_updateFIFObuffer(reference_frames, nRefFrames, reconstructed_frame)
% Update the reference frames buffer in FIFO manner
if length(reference_frames) >= nRefFrames
    reference_frames(1:end-1) = reference_frames(2:end);  % Shift frames left
end
reference_frames{min(length(reference_frames) + 1, nRefFrames)} = reconstructed_frame;