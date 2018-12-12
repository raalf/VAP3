function [] = fcnXMLPANEL(x,y,z,chord,twist,airfoil,filename)
% This function creates a panel for the .vap input file
%
% Inputs:
%   x,y,z - vector of x,y,z coordinates respectively (m)
%   chord - vector of the chord distribution (m)
%   twist - vector of the twist distribtuion (deg)
%   airfoil - string of the desired panel airfoil
%   filename (optional) - add filename to save filename.vap
%
% Outputs:
%   filename.vap file
%   Command window printout of .vap file
%
% 
% Sample function call:
% x = [1 2 3 4 5];
% y = [5 4 3 2 1];
% z = [1 3 1 2 1];
% chord = [0.1 0.1 0.1 0.1 0.5];
% twist = [5 5 5 5 9];
% fcnXMLPANEL(x,y,z,chord,twist,'NACA0012')



% Check to make sure x,y,z,chord,twist are all the same length
if length(x) ~= length(y) || length(x) ~= length(z) || ...
        length(x) ~= length(chord) || length(x) ~= length(twist)
    fprintf('Uh Oh! Your x,y,z,chord,twist inputs dont match dimensions\n')
    return
end

% Create Panel node (main node)
docNode = com.mathworks.xml.XMLUtils.createDocument('panel');
docRootNode = docNode.getDocumentElement;

% Create airfoil, symmetry and N nodes and write their values
thisElement = docNode.createElement('airfoil'); 
thisElement.appendChild(docNode.createTextNode(num2str(airfoil)));
docRootNode.appendChild(thisElement);

thisElement = docNode.createElement('symmetry'); 
thisElement.appendChild(docNode.createTextNode('1'));
docRootNode.appendChild(thisElement);

thisElement = docNode.createElement('N'); 
thisElement.appendChild(docNode.createTextNode('1'));
docRootNode.appendChild(thisElement);

for i=1:length(x)
    % Create section node
    section_node = docNode.createElement('section');
    docNode.getDocumentElement.appendChild(section_node);
    
    x_node = docNode.createElement('x'); % Create x element
    x_node.appendChild(docNode.createTextNode(num2str(x(i)))); %Add x value
    docRootNode.appendChild(x_node); % Add x to xml
    section_node.appendChild(x_node); % Add x under the section

    y_node = docNode.createElement('y');
    y_node.appendChild(docNode.createTextNode(num2str(y(i))));
    docRootNode.appendChild(y_node);
    section_node.appendChild(y_node);

    z_node = docNode.createElement('z');
    z_node.appendChild(docNode.createTextNode(num2str(z(i))));
    docRootNode.appendChild(z_node);
    section_node.appendChild(z_node);

    chord_node = docNode.createElement('chord');
    chord_node.appendChild(docNode.createTextNode(num2str(chord(i))));
    docRootNode.appendChild(chord_node);
    section_node.appendChild(chord_node);

    chord_node = docNode.createElement('twist');
    chord_node.appendChild(docNode.createTextNode(num2str(twist(i))));
    docRootNode.appendChild(chord_node);
    section_node.appendChild(chord_node);
end

% Save to filename if inputed to function otherwise just print to command
% window
if exist('filename','var')
    xmlFileName = [filename,'.vap'];
    xmlwrite(xmlFileName,docNode);
    type(xmlFileName);
else
    xmlFileName = [tempname,'.vap'];
    xmlwrite(xmlFileName,docNode);
    type(xmlFileName);
    delete(xmlFileName)
end
