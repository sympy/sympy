<?xml version="1.0" encoding="utf-8"?>
<!-- ********************************************************** -->
<!-- XSL Transform of OpenMath to Simplified Content MathML     -->
<!--                                                            -->
<!-- Author: Clare M. So <clare@scl.csd.uwo.ca>                 -->
<!--                                                            -->
<!-- September 2002                                             -->
<!-- ********************************************************** -->


<!-- Special MathML Entities -->

<!DOCTYPE stylesheet [
<!ENTITY pi "&#x003C0;">
<!ENTITY e "&#x02147E;">
<!ENTITY ee "&#x02147E;">
<!ENTITY ExponentialE "&#x02147E;">
<!ENTITY ImaginaryI "&#x02148;">
<!ENTITY ii "&#x02148;">
<!ENTITY gamma "&#x003B3;">
<!ENTITY infin "&#x0221E;">
<!ENTITY infty "&#x0221E;">
<!ENTITY true "&#xF0002;">
<!ENTITY false "&#xF0003;">
<!ENTITY NotANumber "&#xF0001;">
<!ENTITY NaN "&#xF0001;">
]>


<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns="http://www.w3.org/1998/Math/MathML"
  xmlns:om="http://www.openmath.org/OpenMath"
  exclude-result-prefixes="om"
  version="1.0">

  <xsl:output method="xml" indent="yes"/>
  
  <xsl:strip-space elements="*"/>


  <!-- OMOBJ (D. Carlisle) -->
  <xsl:template match="om:OMOBJ">
    <math>
      <xsl:apply-templates/>
    </math>
  </xsl:template>


  <!-- OMI (D. Carlisle) -->
  <xsl:template match="om:OMI">
    <cn type="integer">
      <xsl:variable name="x" select="."/>
      <xsl:choose>
        <xsl:when test="contains($x,'x')">
          <xsl:attribute name="base">16</xsl:attribute>
          <xsl:value-of select="normalize-space(concat(substring-before($x,'x'),substring-after($x,'x')))"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="normalize-space($x)"/> <!-- default is decimal -->
        </xsl:otherwise>
      </xsl:choose>
    </cn>
  </xsl:template>


  <!-- OMV (D. Carlisle) -->
  <xsl:template match="om:OMV">
    <ci>
      <xsl:value-of select="normalize-space(@name)"/>
    </ci>
  </xsl:template>


  <!-- OMF (decimal) -->  
  <xsl:template match="om:OMF[@dec]">
    <cn>
      <xsl:if test="contains(@dec,'e')">
        <xsl:attribute name="type">e-natation</xsl:attribute>
      </xsl:if>
      <xsl:value-of select="normalize-space(@dec)"/>
    </cn>
  </xsl:template>
  

  <!-- OMF (hex)-->
  <xsl:template match="om:OMF[@hex]">
    <cn base="16">
      <xsl:if test="contains(@hex,'e')">
        <xsl:attribute name="type">e-natation</xsl:attribute>  
      </xsl:if>
      <xsl:value-of select="normalize-space(@hex)"/>
    </cn>
  </xsl:template>


  <!-- Change OMA to apply (D. Carlisle) -->
  <xsl:template match="om:OMA">
    <apply>
      <xsl:apply-templates/>
    </apply>
  </xsl:template>

  
  <!-- OMB -->
  <xsl:template match="om:OMB">
    <mtext definitionURL="http://www.openmath.org/objects#OMB">
      <xsl:value-of select="."/>
    </mtext>
  </xsl:template>


  <!-- OME -->
  <xsl:template match="om:OME"/>
    

  <!-- OMSTR (D. Carlisle) -->
  <xsl:template match="om:OMSTR">
    <mtext>
      <xsl:value-of select="."/>
    </mtext>
  </xsl:template>
  
  
  <!-- OMATTR (no direct equivalent in MathML) -->
  <xsl:template match="om:OMATTR[om:OMATP[om:OMS[@cd='altenc']]]">
    <semantics>
      <xsl:apply-templates/>
    </semantics>
  </xsl:template>

  <!-- All mathmltypes -->
  <xsl:template match="om:OMATTR[om:OMATP[om:OMS[@cd='mathmltypes']]]">
    <ci>
      <xsl:attribute name="type">
        <xsl:value-of select="normalize-space(translate(substring-before(om:OMATP/om:OMS[2]/@name,'_type'),'_','-'))"/>
      </xsl:attribute>
      <xsl:value-of select="normalize-space(*[2]/@name)"/>
    </ci>
  </xsl:template>


  <!-- OMATP -->
  <xsl:template match="om:OMATP">
    <xsl:apply-templates/>
  </xsl:template>


  <!-- OMBIND -->
  <xsl:template match="om:OMBIND">
    <apply>
      <xsl:apply-templates/>
    </apply>
  </xsl:template>


  <!-- OMBVAR -->
  <xsl:template match="om:OMBVAR">
    <xsl:for-each select="om:OMV">
      <bvar><xsl:apply-templates select="."/></bvar>
    </xsl:for-each>
  </xsl:template>


  <!-- all OMS's -->
  <xsl:template match="om:OMS">
    <csymbol>
      <xsl:attribute name="definitionURL">
        <xsl:value-of select="concat(concat(concat('http://www.openmath.org/cd/',@cd),'#'),@name)"/>
      </xsl:attribute>
    </csymbol>
  </xsl:template>

</xsl:stylesheet>






