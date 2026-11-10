function [] = plt(arg)

arg=squeeze(arg);
% s = size(arg);
% if s(1)==1
%     arg=arg.';
% end


figure();
switch ndims(arg)-sum(ndims(arg)==1)
    case 1
        if isreal(arg)
            plot(arg);legend('real(arg)');
        elseif isreal(imag(arg))
            plot(imag(arg));legend('imag(arg)');
        else
            plot(abs(arg));legend('abs(arg)');
        end
    case 2
        if isreal(arg)
            imagesc(arg);disp('real(arg)');
        elseif isreal(imag(arg))
            imagesc(imag(arg));disp('imag(arg)');
        else
            imagesc(abs(arg));disp('abs(arg)');
        end
    % case 3
    %     disp('dim too high')
end

end

