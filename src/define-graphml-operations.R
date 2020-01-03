# TODO: Add comment
# 
# Author: nejat
###############################################################################


library(package="igraph")
library(package="XML")


########################################################################
# It removes the XML nodes (starting with "<key") which has an "id" attribute mathced with the input parameter "id".
#
# doc: the content of the XML document. It is obtained through "xmlParse(file.path)"
# id: the attribute called "id", which will be looked for inside a XML node called "key"
# 
# EXAMPLE: If the XML file contains only 2 "key" nodes, like below, and if id is matched with one of "e_weight" and "e_sign",
#   it will remove the corresponding "key" node.
# <key id="e_weight" for="edge" attr.name="weight" attr.type="double"/>
# <key id="e_sign" for="edge" attr.name="sign" attr.type="double"/>
#
########################################################################
removeAttrDefinitionNode = function(doc, id){
	
	keys = getNodeSet(doc, "//*[name()='key']")
	nodeToBeRemoved = NA
	
	for(key in keys){
		if(xmlAttrs(key)["id"] == id){
			nodeToBeRemoved = key
			break
		}
	}
	
	if(!is.na(nodeToBeRemoved)){
		removeNodes(nodeToBeRemoved)
	} else{
		print("error: remove nodes")
	}
}



########################################################################
# It adds a new XML tag, i.e. attribute definition node, just before the last XML node corresponding to the graph content,
#   then returns the content of the updated graphml file
#
# doc: the content of the XML document. It is obtained through "xmlParse(file.path)"
# id: a new attribute to be created in which the name is "id"
# for.name: a new attribute to be created in which the name is "for"
# attr.name: a new attribute to be created in which the name is "attr.name"
# attr.type: a new attribute to be created in which the name is "attr.type"
#
#
# EXEMPLE:
# addAttrDefinitionNode(doc, id="e_sign", for.name="edge", attr.name="sign", attr.type="double")
# The result: <key id="e_sign" for="edge" attr.name="sign" attr.type="double"/>
########################################################################
addAttrDefinitionNode = function(doc, id, for.name, attr.name, attr.type){
	
	x = xmlRoot(doc, skip = FALSE)
	s = newXMLNode("key", parent=x)
	addAttributes(s, "id"=id)
	addAttributes(s, "for"=for.name)
	addAttributes(s, "attr.name"=attr.name)
	addAttributes(s, "attr.type"=attr.type)
	
	# Before adding a new definition node, the last child of the root node had the graph content
	#  (nodes def, edges def) which has many many lines.
	# And adding a new definition node with "newXMLNode()" puts the XML tage in the end
	# To overcome this situation, i.e. to keep this graph content at the last XML
	#  do the following:
	n=xmlSize(x)
	last.before.node = x[[n-1]] # originally the graph content
	x[[n-1]] = x[[n]] # x[[n]] is new XML node
	x[[n]] = last.before.node
	
}


########################################################################
# It adds partition information conveyed by the argument "mems". To store this information, an attribute name should be given.
#
# file.path: file path for the XML file
# attr.name: a value to be used for the attribute name. Example: "solution23"
# mems: partition information, associated with "attr.name". It is a vector
#
########################################################################
addPartitionInfoIntoNodes = function(file.path, attr.name, mems)
{
	doc <- xmlParse(file.path)
	
	id = paste("v", attr.name, sep="_")
	addAttrDefinitionNode(
			doc, id=id, for.name="node", attr.name=attr.name, attr.type="double"
	)
	
	nodes = getNodeSet(doc, "//*[name()='node']")
	for(i in 1:length(nodes)){
		node = nodes[i]
		node.mem.info = mems[i]
		partition = newXMLNode("data", parent=node)
		addAttributes(partition, "key"=id)
		xmlValue(partition) = node.mem.info
	}
	
	return(doc)
}



########################################################################
# It adds the sign attributes into the edges. It takes the file path as input argument, 
#    then calls the method "addSignAttrIntoDocument"
#
########################################################################
addSignAttrIntoGraphmlFile = function(file.path){
    doc <- xmlParse(file.path)
    addSignAttrIntoDocument(doc)
}

########################################################################
# It adds the sign attributes into the edges. It takes the XML content as input argument.
#
# When a sign attribute is adeded into an edge, its content becomes:
# <edge source="n3" target="n4">
#   <data key="e_weight">1</data>
#    <data key="e_sign">1</data>
# </edge>
#
# doc: the content of the XML document. It is obtained through "xmlParse(file.path)"
#   
########################################################################
addSignAttrIntoDocument = function(doc){

	addAttrDefinitionNode(
			doc, id="e_sign", for.name="edge", attr.name="sign", attr.type="double"
	)
	
	edges = getNodeSet(doc, "//*[name()='edge']")
	for(edge in edges){
		e = xmlChildren(edge)
		e.data = e[[1]] # first child
#		print(e.data)
		attr.name = xmlAttrs(e.data)["key"]
		if(attr.name == "e_weight" || attr.name == "weight" ){
			val = as.numeric( xmlValue(e.data) )
			xmlValue(e.data) = abs(val)
			s = newXMLNode("data", parent=edge)
			addAttributes(s, "key"="e_sign")
			if(val<0){	
				xmlValue(s) = -1
			} else{
				xmlValue(s) = 1	
			}
		}
	}
	
	return(doc)
}


########################################################################
# It removes the sign attributes from the edges: For remainder, when a sign attribute is adeded into an edge, its content becomes:
# <edge source="n3" target="n4">
#   <data key="e_weight">1</data>
#    <data key="e_sign">1</data>
# </edge>
#
# file.path: file path for the XML file
#       
########################################################################
removeSignAttrFromGraphmlFile = function(file.path){
	doc <- xmlParse(file.path)
	
	removeAttrDefinitionNode(doc, "e_sign")
	
	edges = getNodeSet(doc, "//*[name()='edge']")
	for(edge in edges){
		e = xmlChildren(edge)
		e.weight = e[[1]] # first child
		e.sign = e[[2]] # second child
		val = as.numeric( xmlValue(e.weight) )
		if(as.numeric(xmlValue(e.sign)) == -1){
			xmlValue(e.weight) = val * -1
		}
		
		removeNodes(e.sign)
		
	}	
	
	
	return(doc)
}





# ==============================================================================

# EXAMPLE.GRAPHML.CONTENT = '<?xml version="1.0" encoding="UTF-8"?>
# <graphml xmlns="http://graphml.graphdrawing.org/xmlns"
#      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
#      xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
#     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd"> 
#   <key id="e_weight" for="edge" attr.name="weight" attr.type="double"/>
#   <graph id="G" edgedefault="undirected">
#     <node id="n0">
#     </node>
#     <node id="n1">
#     </node>
#     <node id="n2">
#     </node>
#     <node id="n3">
#     </node>
#     <node id="n4">
#     </node>
#     <edge source="n0" target="n1">
#       <data key="e_weight">1</data>
#     </edge>
#     <edge source="n0" target="n2">
#       <data key="e_weight">1</data>
#     </edge>
#     <edge source="n0" target="n3">
#       <data key="e_weight">1</data>
#     </edge>
#     <edge source="n0" target="n4">
#       <data key="e_weight">-1</data>
#     </edge>
#     <edge source="n1" target="n2">
#       <data key="e_weight">1</data>
#     </edge>
#     <edge source="n1" target="n3">
#       <data key="e_weight">1</data>
#     </edge>
#     <edge source="n1" target="n4">
#       <data key="e_weight">1</data>
#     </edge>
#     <edge source="n2" target="n3">
#       <data key="e_weight">1</data>
#     </edge>
#     <edge source="n2" target="n4">
#       <data key="e_weight">-1</data>
#     </edge>
#     <edge source="n3" target="n4">
#       <data key="e_weight">1</data>
#     </edge>
#   </graph>
# </graphml>
# '
# 
# # ==========
# # Scenario 1 
# # ==========
# 
# file.path = "temp-graphml-ops.graphml"
# write(x=EXAMPLE.GRAPHML.CONTENT, file=file.path)
# 
# doc = addSignAttrIntoGraphmlFile(file.path)
# print(doc)
# 
# 
# # ------
# 
# # ==========
# # Scenario 2
# # ==========
# 
# file.path = "temp-graphml-ops.graphml"
# write(x=EXAMPLE.GRAPHML.CONTENT, file=file.path)
# 
# nodes = getNodeSet(doc, "//*[name()='node']")
# 
# addAttrDefinitionNode(
#     doc, id="e_sign", for.name="edge", attr.name="sign", attr.type="double"
# )
# 
# keys = getNodeSet(doc, "//*[name()='key']")
# for(key in keys)
#     print(xmlAttrs(key)["id"])
# 
# print(doc)
# removeAttrDefinitionNode(doc, "e_sign")
# print(doc)

