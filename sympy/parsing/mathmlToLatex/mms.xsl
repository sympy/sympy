<!-- Took reference from github repository : https://github.com/oerpub/mathconverter/tree/master/xsl_yarosh -->

<?xml version='1.0' encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
		xmlns:m="http://www.w3.org/1998/Math/MathML"
                version='1.0'>
                
<xsl:output method="text" indent="no" encoding="UTF-8"/>

<!-- ====================================================================== -->
<!-- $Id: mmltex.xsl,v 1.7 2003/06/10 12:24:04 shade33 Exp $
     This file is part of the XSLT MathML Library distribution.
     See ./README or http://www.raleigh.ru/MathML/mmltex for
     copyright and other information                                        -->
<!-- ====================================================================== -->

<xsl:include href="tokens.xsl"/>
<xsl:include href="glayout.xsl"/>
<xsl:include href="scripts.xsl"/>
<xsl:include href="tables.xsl"/>
<xsl:include href="entities.xsl"/>
<xsl:include href="cmarkup.xsl"/>

<xsl:strip-space elements="m:*"/>

<xsl:template match="m:math[not(@mode) or @mode='inline'][not(@display)] | m:math[@display='inline']">
	<xsl:text>&#x00024; </xsl:text>
	<xsl:apply-templates/>
	<xsl:text>&#x00024;</xsl:text>
</xsl:template>

<xsl:template match="m:math[@display='block'] | m:math[@mode='display'][not(@display)]">
	<xsl:text>&#xA;\[&#xA;&#x9;</xsl:text>
	<xsl:apply-templates/>
	<xsl:text>&#xA;\]</xsl:text>
</xsl:template>

</xsl:stylesheet>