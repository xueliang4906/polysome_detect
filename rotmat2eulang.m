function eul = rotmat2eulang(R,order)
% ROTMAT2EULANG Calculates fixed-euler angles from rotation matrix
% 20201030 download from https://github.com/wspr/matlab-euler-angles/blob/master/rotmat2eulang.m
% * `R`: Rotation matrix (3x3)
% * `order`: Order of Euler angles to return; any of:
%              {'XYZ','XZY','YXZ','YZX','ZYX','ZXY'}
% * `EUL`: Euler angles in specified order
%

switch order

  case 'XYZ'
    eul(1) = atan2(-R(2,3),R(3,3)).*180./pi;
    eul(2) = asin(R(1,3)).*180./pi;
    eul(3) = atan2(-R(1,2),R(1,1)).*180./pi;
    
  case 'XZY'
    eul(1) = atan2(R(3,2),R(2,2)).*180./pi;
    eul(2) = -asin(R(1,2)).*180./pi;
    eul(3) = atan2(R(1,3),R(1,1)).*180./pi;
    
  case 'YXZ'
    eul(1) = atan2(R(1,3),R(3,3)).*180./pi;
    eul(2) = -asin(R(2,3)).*180./pi;
    eul(3) = atan2(R(2,1),R(2,2)).*180./pi;
    
  case 'YZX'
    eul(1) = atan2(-R(3,1),R(1,1)).*180./pi;
    eul(2) = asin(R(2,1)).*180./pi;
    eul(3) = atan2(-R(2,3),R(2,2)).*180./pi;
    
  case 'ZXY'
    eul(1) = atan2(-R(1,2),R(2,2)).*180./pi;
    eul(2) = asin(R(3,2)).*180./pi;
    eul(3) = atan2(-R(3,1),R(3,3)).*180./pi;

  case 'ZYX'
    eul(1) = atan2(R(2,1),R(1,1)).*180./pi;
    eul(2) = -asin(R(3,1)).*180./pi;
    eul(3) = atan2(R(3,2),R(3,3)).*180./pi;
    
  otherwise
    error('Don''t know specified rotation order %s', order)

end

end

% Licence included in README.