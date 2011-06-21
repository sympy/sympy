
    

  

<!DOCTYPE html>
<html>
  <head>
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="chrome=1">
    <script type="text/javascript">var NREUMQ=[];NREUMQ.push(["mark","firstbyte",new Date().getTime()]);</script>
        <title>sympy/statistics/testcase.py at fit_functions_01 from andymd/sympy - GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />

    <link href="https://a248.e.akamai.net/assets.github.com/3caf378f22640d86d06e445746d14a0fb0742a68/stylesheets/bundle_github.css" media="screen" rel="stylesheet" type="text/css" />
    

    <script type="text/javascript">
      if (typeof console == "undefined" || typeof console.log == "undefined")
        console = { log: function() {} }
    </script>
    <script type="text/javascript" charset="utf-8">
      var GitHub = {
        assetHost: 'https://a248.e.akamai.net/assets.github.com'
      }
      var github_user = 'andymd'
      
    </script>
    <script src="https://a248.e.akamai.net/assets.github.com/javascripts/jquery/jquery-1.6.1.min.js" type="text/javascript"></script>
    <script src="https://a248.e.akamai.net/assets.github.com/028b735b7be9bb4df802d5efba13d1360dfe7810/javascripts/bundle_github.js" type="text/javascript"></script>


    
    <script type="text/javascript" charset="utf-8">
      GitHub.spy({
        repo: "andymd/sympy"
      })
    </script>

    
  <link rel='canonical' href='/andymd/sympy/blob/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf/sympy/statistics/testcase.py'>

  <link href="https://github.com/andymd/sympy/commits/fit_functions_01.atom" rel="alternate" title="Recent Commits to sympy:fit_functions_01" type="application/atom+xml" />

    <META NAME="ROBOTS" CONTENT="NOINDEX, FOLLOW">    <meta name="description" content="sympy - A computer algebra system written in pure Python" />
    <script type="text/javascript">
      GitHub.nameWithOwner = GitHub.nameWithOwner || "andymd/sympy";
      GitHub.currentRef = 'fit_functions_01';
      GitHub.commitSHA = "9475e4d85aa1c917dbf31bc8a38abad52b5ecadf";
      GitHub.currentPath = 'sympy/statistics/testcase.py';
      GitHub.masterBranch = "master";

      
    </script>
  

        <script type="text/javascript">
      var _gaq = _gaq || [];
      _gaq.push(['_setAccount', 'UA-3769691-2']);
      _gaq.push(['_setDomainName', 'none']);
      _gaq.push(['_trackPageview']);
      _gaq.push(['_trackPageLoadTime']);
      (function() {
        var ga = document.createElement('script');
        ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
        ga.setAttribute('async', 'true');
        document.documentElement.firstChild.appendChild(ga);
      })();
    </script>

    
  </head>

  

  <body class="logged_in page-blob  linux env-production">
    

    

    

    <div class="subnavd" id="main">
      <div id="header" class="true">
        
          <a class="logo boring" href="https://github.com/">
            
            <img alt="github" class="default" height="45" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov5.png" />
            <!--[if (gt IE 8)|!(IE)]><!-->
            <img alt="github" class="hover" height="45" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov5-hover.png" />
            <!--<![endif]-->
          </a>
        
        
          





  
    <div class="userbox">
      <div class="avatarname">
        <a href="https://github.com/andymd"><img src="https://secure.gravatar.com/avatar/15e0232481d6fabaa32a5421a0224398?s=140&d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-140.png" alt="" width="20" height="20"  /></a>
        <a href="https://github.com/andymd" class="name">andymd</a>

        
        
          <a href="https://github.com/inbox/notifications" class="unread_count notifications_count new tooltipped downwards js-notification-count" title="Unread Notifications">3</a>
        
      </div>
      <ul class="usernav">
        <li><a href="https://github.com/">Dashboard</a></li>
        <li>
          
          <a href="https://github.com/inbox">Inbox</a>
          <a href="https://github.com/inbox" class="unread_count ">0</a>
        </li>
        <li><a href="https://github.com/account">Account Settings</a></li>
        <li><a href="/logout">Log Out</a></li>
      </ul>
    </div><!-- /.userbox -->
  


        
        <div class="topsearch">
  
    <form action="/search" id="top_search_form" method="get">
      <a href="/search" class="advanced-search tooltipped downwards" title="Advanced Search">Advanced Search</a>
      <div class="search placeholder-field js-placeholder-field">
        <label class="placeholder" for="global-search-field">Search…</label>
        <input type="text" class="search my_repos_autocompleter" id="global-search-field" name="q" results="5" /> <input type="submit" value="Search" class="button" />
      </div>
      <input type="hidden" name="type" value="Everything" />
      <input type="hidden" name="repo" value="" />
      <input type="hidden" name="langOverride" value="" />
      <input type="hidden" name="start_value" value="1" />
    </form>
    <ul class="nav">
      <li><a href="/explore">Explore GitHub</a></li>
      <li><a href="https://gist.github.com">Gist</a></li>
      
      <li><a href="/blog">Blog</a></li>
      
      <li><a href="http://help.github.com">Help</a></li>
    </ul>
  
</div>

      </div>

      
      
        
    <div class="site">
      <div class="pagehead repohead vis-public fork   instapaper_ignore readability-menu">

      

      <div class="title-actions-bar">
        <h1>
          <a href="/andymd">andymd</a> / <strong><a href="/andymd/sympy">sympy</a></strong>
          
            <span class="fork-flag">
              
              <span class="text">forked from <a href="/sympy/sympy">sympy/sympy</a></span>
            </span>
          
          
        </h1>

        
    <ul class="actions">
      

      
        <li class="for-owner" style="display:none"><a href="/andymd/sympy/admin" class="minibutton btn-admin "><span><span class="icon"></span>Admin</span></a></li>
        <li>
          <a href="/andymd/sympy/toggle_watch" class="minibutton btn-watch " id="watch_button" onclick="var f = document.createElement('form'); f.style.display = 'none'; this.parentNode.appendChild(f); f.method = 'POST'; f.action = this.href;var s = document.createElement('input'); s.setAttribute('type', 'hidden'); s.setAttribute('name', 'authenticity_token'); s.setAttribute('value', 'a58c21d49dbae0a97a9c153e4e56a698518b585c'); f.appendChild(s);f.submit();return false;" style="display:none"><span><span class="icon"></span>Watch</span></a>
          <a href="/andymd/sympy/toggle_watch" class="minibutton btn-watch " id="unwatch_button" onclick="var f = document.createElement('form'); f.style.display = 'none'; this.parentNode.appendChild(f); f.method = 'POST'; f.action = this.href;var s = document.createElement('input'); s.setAttribute('type', 'hidden'); s.setAttribute('name', 'authenticity_token'); s.setAttribute('value', 'a58c21d49dbae0a97a9c153e4e56a698518b585c'); f.appendChild(s);f.submit();return false;" style="display:none"><span><span class="icon"></span>Unwatch</span></a>
        </li>
        
          
            <li class="for-notforked" style="display:none"><a href="/andymd/sympy/fork" class="minibutton btn-fork " id="fork_button" onclick="var f = document.createElement('form'); f.style.display = 'none'; this.parentNode.appendChild(f); f.method = 'POST'; f.action = this.href;var s = document.createElement('input'); s.setAttribute('type', 'hidden'); s.setAttribute('name', 'authenticity_token'); s.setAttribute('value', 'a58c21d49dbae0a97a9c153e4e56a698518b585c'); f.appendChild(s);f.submit();return false;"><span><span class="icon"></span>Fork</span></a></li>
            <li class="for-hasfork" style="display:none"><a href="#" class="minibutton btn-fork " id="your_fork_button"><span><span class="icon"></span>Your Fork</span></a></li>
          

          <li id='pull_request_item' class='nspr' style='display:none'><a href="/andymd/sympy/pull/new/fit_functions_01" class="minibutton btn-pull-request "><span><span class="icon"></span>Pull Request</span></a></li>
        
      
      
      <li class="repostats">
        <ul class="repo-stats">
          <li class="watchers"><a href="/andymd/sympy/watchers" title="Watchers" class="tooltipped downwards">1</a></li>
          <li class="forks"><a href="/andymd/sympy/network" title="Forks" class="tooltipped downwards">103</a></li>
        </ul>
      </li>
    </ul>

      </div>

        
  <ul class="tabs">
    <li><a href="/andymd/sympy/tree/fit_functions_01" class="selected" highlight="repo_source">Source</a></li>
    <li><a href="/andymd/sympy/commits/fit_functions_01" highlight="repo_commits">Commits</a></li>
    <li><a href="/andymd/sympy/network" highlight="repo_network">Network</a></li>
    <li><a href="/andymd/sympy/pulls" highlight="repo_pulls">Pull Requests (0)</a></li>

    
      <li><a href="/andymd/sympy/forkqueue" highlight="repo_fork_queue">Fork Queue</a></li>
    

    

                <li><a href="/andymd/sympy/wiki" highlight="repo_wiki">Wiki (0)</a></li>
        
    <li><a href="/andymd/sympy/graphs" highlight="repo_graphs">Graphs</a></li>

    <li class="contextswitch nochoices">
      <span class="toggle leftwards" >
        <em>Branch:</em>
        <code>fit_functions_01</code>
      </span>
    </li>
  </ul>

  <div style="display:none" id="pl-description"><p><em class="placeholder">click here to add a description</em></p></div>
  <div style="display:none" id="pl-homepage"><p><em class="placeholder">click here to add a homepage</em></p></div>

  <div class="subnav-bar">
  
  <ul>
    <li>
      <a href="/andymd/sympy/branches" class="dropdown">Switch Branches (4)</a>
      <ul>
        
          
          
            <li><a href="/andymd/sympy/blob/0.7.0/sympy/statistics/testcase.py" data-name="0.7.0">0.7.0</a></li>
          
        
          
            <li><strong>fit_functions_01 &#x2713;</strong></li>
            
          
          
            <li><a href="/andymd/sympy/blob/gh-pages/sympy/statistics/testcase.py" data-name="gh-pages">gh-pages</a></li>
          
        
          
          
            <li><a href="/andymd/sympy/blob/master/sympy/statistics/testcase.py" data-name="master">master</a></li>
          
        
      </ul>
    </li>
    <li>
      <a href="#" class="dropdown ">Switch Tags (34)</a>
              <ul>
                      
              <li><a href="/andymd/sympy/blob/sympy-0.7.0.rc1/sympy/statistics/testcase.py" data-name="sympy-0.7.0.rc1">sympy-0.7.0.rc1</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.7/sympy/statistics/testcase.py" data-name="sympy-0.6.7">sympy-0.6.7</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.6.rc1/sympy/statistics/testcase.py" data-name="sympy-0.6.6.rc1">sympy-0.6.6.rc1</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.6/sympy/statistics/testcase.py" data-name="sympy-0.6.6">sympy-0.6.6</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.5.rc2/sympy/statistics/testcase.py" data-name="sympy-0.6.5.rc2">sympy-0.6.5.rc2</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.5.rc1/sympy/statistics/testcase.py" data-name="sympy-0.6.5.rc1">sympy-0.6.5.rc1</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.5.beta3/sympy/statistics/testcase.py" data-name="sympy-0.6.5.beta3">sympy-0.6.5.beta3</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.5.beta2/sympy/statistics/testcase.py" data-name="sympy-0.6.5.beta2">sympy-0.6.5.beta2</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.5-beta1/sympy/statistics/testcase.py" data-name="sympy-0.6.5-beta1">sympy-0.6.5-beta1</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.5/sympy/statistics/testcase.py" data-name="sympy-0.6.5">sympy-0.6.5</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.4.beta3/sympy/statistics/testcase.py" data-name="sympy-0.6.4.beta3">sympy-0.6.4.beta3</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.4.beta2/sympy/statistics/testcase.py" data-name="sympy-0.6.4.beta2">sympy-0.6.4.beta2</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.4.beta1/sympy/statistics/testcase.py" data-name="sympy-0.6.4.beta1">sympy-0.6.4.beta1</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.4/sympy/statistics/testcase.py" data-name="sympy-0.6.4">sympy-0.6.4</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.3.beta2/sympy/statistics/testcase.py" data-name="sympy-0.6.3.beta2">sympy-0.6.3.beta2</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.3.beta1/sympy/statistics/testcase.py" data-name="sympy-0.6.3.beta1">sympy-0.6.3.beta1</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.3/sympy/statistics/testcase.py" data-name="sympy-0.6.3">sympy-0.6.3</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.2/sympy/statistics/testcase.py" data-name="sympy-0.6.2">sympy-0.6.2</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.1/sympy/statistics/testcase.py" data-name="sympy-0.6.1">sympy-0.6.1</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.6.0/sympy/statistics/testcase.py" data-name="sympy-0.6.0">sympy-0.6.0</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.15/sympy/statistics/testcase.py" data-name="sympy-0.5.15">sympy-0.5.15</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.14/sympy/statistics/testcase.py" data-name="sympy-0.5.14">sympy-0.5.14</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.13/sympy/statistics/testcase.py" data-name="sympy-0.5.13">sympy-0.5.13</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.12/sympy/statistics/testcase.py" data-name="sympy-0.5.12">sympy-0.5.12</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.11/sympy/statistics/testcase.py" data-name="sympy-0.5.11">sympy-0.5.11</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.10/sympy/statistics/testcase.py" data-name="sympy-0.5.10">sympy-0.5.10</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.9/sympy/statistics/testcase.py" data-name="sympy-0.5.9">sympy-0.5.9</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.8/sympy/statistics/testcase.py" data-name="sympy-0.5.8">sympy-0.5.8</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.7/sympy/statistics/testcase.py" data-name="sympy-0.5.7">sympy-0.5.7</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.6/sympy/statistics/testcase.py" data-name="sympy-0.5.6">sympy-0.5.6</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.5/sympy/statistics/testcase.py" data-name="sympy-0.5.5">sympy-0.5.5</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.4/sympy/statistics/testcase.py" data-name="sympy-0.5.4">sympy-0.5.4</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.3/sympy/statistics/testcase.py" data-name="sympy-0.5.3">sympy-0.5.3</a></li>
            
                      
              <li><a href="/andymd/sympy/blob/sympy-0.5.0/sympy/statistics/testcase.py" data-name="sympy-0.5.0">sympy-0.5.0</a></li>
            
                  </ul>
      
    </li>
    <li>
    
    <a href="/andymd/sympy/branches/fit_functions_01" class="manage">Branch List</a>
    
    </li>
  </ul>
</div>

  
  
  
  
  
  



        
    <div id="repo_details" class="metabox clearfix">
      <div id="repo_details_loader" class="metabox-loader" style="display:none">Sending Request&hellip;</div>

        <a href="/andymd/sympy/downloads" class="download-source" id="download_button" title="Download source, tagged packages and binaries."><span class="icon"></span>Downloads</a>

      <div id="repository_desc_wrapper">
      <div id="repository_description" rel="repository_description_edit">
        
          <p>A computer algebra system written in pure Python
            <span id="read_more" style="display:none">&mdash; <a href="#readme">Read more</a></span>
          </p>
        
      </div>

      <div id="repository_description_edit" style="display:none;" class="inline-edit">
        <form action="/andymd/sympy/admin/update" method="post"><div style="margin:0;padding:0"><input name="authenticity_token" type="hidden" value="a58c21d49dbae0a97a9c153e4e56a698518b585c" /></div>
          <input type="hidden" name="field" value="repository_description">
          <input type="text" class="textfield" name="value" value="A computer algebra system written in pure Python">
          <div class="form-actions">
            <button class="minibutton"><span>Save</span></button> &nbsp; <a href="#" class="cancel">Cancel</a>
          </div>
        </form>
      </div>

      
      <div class="repository-homepage" id="repository_homepage" rel="repository_homepage_edit">
        <p><a href="http://sympy.org/" rel="nofollow">http://sympy.org/</a></p>
      </div>

      <div id="repository_homepage_edit" style="display:none;" class="inline-edit">
        <form action="/andymd/sympy/admin/update" method="post"><div style="margin:0;padding:0"><input name="authenticity_token" type="hidden" value="a58c21d49dbae0a97a9c153e4e56a698518b585c" /></div>
          <input type="hidden" name="field" value="repository_homepage">
          <input type="text" class="textfield" name="value" value="http://sympy.org/">
          <div class="form-actions">
            <button class="minibutton"><span>Save</span></button> &nbsp; <a href="#" class="cancel">Cancel</a>
          </div>
        </form>
      </div>
      </div>
      <div class="rule "></div>
      <div id="url_box" class="url-box">
  

  <ul class="clone-urls">
    
      
        <li id="private_clone_url"><a href="git@github.com:andymd/sympy.git" data-permissions="Read+Write">SSH</a></li>
      
      <li id="http_clone_url"><a href="https://andymd@github.com/andymd/sympy.git" data-permissions="Read+Write">HTTP</a></li>
      <li id="public_clone_url"><a href="git://github.com/andymd/sympy.git" data-permissions="Read-Only">Git Read-Only</a></li>
    
    
  </ul>
  <input type="text" spellcheck="false" id="url_field" class="url-field" />
        <span style="display:none" id="url_box_clippy"></span>
      <span id="clippy_tooltip_url_box_clippy" class="clippy-tooltip tooltipped" title="copy to clipboard">
      <object classid="clsid:d27cdb6e-ae6d-11cf-96b8-444553540000"
              width="14"
              height="14"
              class="clippy"
              id="clippy" >
      <param name="movie" value="https://a248.e.akamai.net/assets.github.com/flash/clippy.swf?v5"/>
      <param name="allowScriptAccess" value="always" />
      <param name="quality" value="high" />
      <param name="scale" value="noscale" />
      <param NAME="FlashVars" value="id=url_box_clippy&amp;copied=&amp;copyto=">
      <param name="bgcolor" value="#FFFFFF">
      <param name="wmode" value="opaque">
      <embed src="https://a248.e.akamai.net/assets.github.com/flash/clippy.swf?v5"
             width="14"
             height="14"
             name="clippy"
             quality="high"
             allowScriptAccess="always"
             type="application/x-shockwave-flash"
             pluginspage="http://www.macromedia.com/go/getflashplayer"
             FlashVars="id=url_box_clippy&amp;copied=&amp;copyto="
             bgcolor="#FFFFFF"
             wmode="opaque"
      />
      </object>
      </span>

  <p id="url_description"><strong>Read+Write</strong> access</p>
</div>

    </div>

    <div class="frame frame-center tree-finder" style="display:none">
      <div class="breadcrumb">
        <b><a href="/andymd/sympy">sympy</a></b> /
        <input class="tree-finder-input" type="text" name="query" autocomplete="off" spellcheck="false">
      </div>

      
        <div class="octotip">
          <p>
            <a href="/andymd/sympy/dismiss-tree-finder-help" class="dismiss js-dismiss-tree-list-help" title="Hide this notice forever">Dismiss</a>
            <strong>Octotip:</strong> You've activated the <em>file finder</em> by pressing <span class="kbd">t</span>
            Start typing to filter the file list. Use <span class="kbd badmono">↑</span> and <span class="kbd badmono">↓</span> to navigate,
            <span class="kbd">enter</span> to view files.
          </p>
        </div>
      

      <table class="tree-browser" cellpadding="0" cellspacing="0">
        <tr class="js-header"><th>&nbsp;</th><th>name</th></tr>
        <tr class="js-no-results no-results" style="display: none">
          <th colspan="2">No matching files</th>
        </tr>
        <tbody class="js-results-list">
        </tbody>
      </table>
    </div>

    <div id="jump-to-line" style="display:none">
      <h2>Jump to Line</h2>
      <form>
        <input class="textfield" type="text">
        <div class="full-button">
          <button type="submit" class="classy">
            <span>Go</span>
          </button>
        </div>
      </form>
    </div>


        

      </div><!-- /.pagehead -->

      

  







<script type="text/javascript">
  GitHub.downloadRepo = '/andymd/sympy/archives/fit_functions_01'
  GitHub.revType = "ref"

  GitHub.repoName = "sympy"
  GitHub.controllerName = "blob"
  GitHub.actionName     = "show"
  GitHub.currentAction  = "blob#show"

  
    GitHub.loggedIn = true
    GitHub.hasWriteAccess = true
    GitHub.hasAdminAccess = true
    GitHub.watchingRepo = true
    GitHub.ignoredRepo = false
    GitHub.hasForkOfRepo = null
    GitHub.hasForked = false
  

  
</script>




<div class="flash-messages"></div>


  <div id="commit">
    <div class="group">
        
  <div class="envelope commit">
    <div class="human">
      
        <div class="message"><pre><a href="/andymd/sympy/commit/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf">added t student distribution</a> </pre></div>
      

      <div class="actor">
        <div class="gravatar">
          
          <img src="https://secure.gravatar.com/avatar/15e0232481d6fabaa32a5421a0224398?s=140&d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-140.png" alt="" width="30" height="30"  />
        </div>
        <div class="name"><a href="/andymd">andymd</a> <span>(author)</span></div>
        <div class="date">
          <time class="js-relative-date" datetime="2011-06-20T08:32:20-07:00" title="2011-06-20 08:32:20">June 20, 2011</time>
        </div>
      </div>

      

    </div>
    <div class="machine">
      <span>c</span>ommit&nbsp;&nbsp;<a href="/andymd/sympy/commit/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf" hotkey="c">9475e4d85aa1c917dbf3</a><br />
      <span>t</span>ree&nbsp;&nbsp;&nbsp;&nbsp;<a href="/andymd/sympy/tree/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf" hotkey="t">6dbd889886b7e15db5a6</a><br />
      
        <span>p</span>arent&nbsp;
        
        <a href="/andymd/sympy/tree/9f41804622c8a012b920cdfe87b2d1c75f31e480" hotkey="p">9f41804622c8a012b920</a>
      

    </div>
  </div>

    </div>
  </div>



  <div id="slider">

  

    <div class="breadcrumb" data-path="sympy/statistics/testcase.py/">
      <b><a href="/andymd/sympy/tree/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf">sympy</a></b> / <a href="/andymd/sympy/tree/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf/sympy">sympy</a> / <a href="/andymd/sympy/tree/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf/sympy/statistics">statistics</a> / testcase.py       <span style="display:none" id="clippy_2069">sympy/statistics/testcase.py</span>
      
      <object classid="clsid:d27cdb6e-ae6d-11cf-96b8-444553540000"
              width="110"
              height="14"
              class="clippy"
              id="clippy" >
      <param name="movie" value="https://a248.e.akamai.net/assets.github.com/flash/clippy.swf?v5"/>
      <param name="allowScriptAccess" value="always" />
      <param name="quality" value="high" />
      <param name="scale" value="noscale" />
      <param NAME="FlashVars" value="id=clippy_2069&amp;copied=copied!&amp;copyto=copy to clipboard">
      <param name="bgcolor" value="#FFFFFF">
      <param name="wmode" value="opaque">
      <embed src="https://a248.e.akamai.net/assets.github.com/flash/clippy.swf?v5"
             width="110"
             height="14"
             name="clippy"
             quality="high"
             allowScriptAccess="always"
             type="application/x-shockwave-flash"
             pluginspage="http://www.macromedia.com/go/getflashplayer"
             FlashVars="id=clippy_2069&amp;copied=copied!&amp;copyto=copy to clipboard"
             bgcolor="#FFFFFF"
             wmode="opaque"
      />
      </object>
      

    </div>

    <div class="frames">
      <div class="frame frame-center" data-path="sympy/statistics/testcase.py/" data-canonical-url="/andymd/sympy/blob/9475e4d85aa1c917dbf31bc8a38abad52b5ecadf/sympy/statistics/testcase.py">
        
          <ul class="big-actions">
            
            <li><a class="file-edit-link minibutton" href="/andymd/sympy/edit/__current_ref__/sympy/statistics/testcase.py"><span>Edit this file</span></a></li>
          </ul>
        

        <div id="files">
          <div class="file">
            <div class="meta">
              <div class="info">
                <span class="icon"><img alt="Txt" height="16" src="https://a248.e.akamai.net/assets.github.com/images/icons/txt.png" width="16" /></span>
                <span class="mode" title="File Mode">100644</span>
                
                  <span>16 lines (12 sloc)</span>
                
                <span>0.338 kb</span>
              </div>
              <ul class="actions">
                <li><a href="/andymd/sympy/raw/fit_functions_01/sympy/statistics/testcase.py" id="raw-url">raw</a></li>
                
                  <li><a href="/andymd/sympy/blame/fit_functions_01/sympy/statistics/testcase.py">blame</a></li>
                
                <li><a href="/andymd/sympy/commits/fit_functions_01/sympy/statistics/testcase.py">history</a></li>
              </ul>
            </div>
            
  <div class="data type-python">
    
      <table cellpadding="0" cellspacing="0">
        <tr>
          <td>
            <pre class="line_numbers"><span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
</pre>
          </td>
          <td width="100%">
            
              
                <div class="highlight"><pre><div class='line' id='LC1'><span class="kn">from</span> <span class="nn">sympy</span> <span class="kn">import</span> <span class="o">*</span></div><div class='line' id='LC2'><span class="kn">from</span> <span class="nn">rseries</span> <span class="kn">import</span> <span class="n">rseries</span></div><div class='line' id='LC3'><span class="kn">from</span> <span class="nn">distributions</span> <span class="kn">import</span> <span class="n">studentt_pdf</span></div><div class='line' id='LC4'><br/></div><div class='line' id='LC5'><span class="n">t</span> <span class="o">=</span> <span class="n">symbols</span><span class="p">(</span><span class="s">&#39;t&#39;</span><span class="p">,</span><span class="n">real</span><span class="o">=</span><span class="bp">True</span><span class="p">)</span></div><div class='line' id='LC6'><span class="n">m</span> <span class="o">=</span> <span class="n">symbols</span><span class="p">(</span><span class="s">&#39;m&#39;</span><span class="p">,</span><span class="n">integer</span><span class="o">=</span><span class="bp">True</span><span class="p">,</span><span class="n">positive</span><span class="o">=</span><span class="bp">True</span><span class="p">)</span></div><div class='line' id='LC7'><br/></div><div class='line' id='LC8'><span class="n">stud</span> <span class="o">=</span> <span class="n">studentt_pdf</span><span class="p">(</span><span class="n">m</span><span class="p">,</span><span class="n">t</span><span class="p">)</span></div><div class='line' id='LC9'><span class="c">#stud = stud.subs(m,10.0)</span></div><div class='line' id='LC10'><br/></div><div class='line' id='LC11'><span class="c">#s,o = rseries( Integral(exp(t),t), t, 0, 50 )</span></div><div class='line' id='LC12'><span class="k">print</span> <span class="s">&quot;_&quot;</span><span class="o">*</span><span class="mi">50</span></div><div class='line' id='LC13'><span class="n">s</span><span class="p">,</span><span class="n">o</span> <span class="o">=</span> <span class="n">rseries</span><span class="p">(</span> <span class="n">Integral</span><span class="p">(</span><span class="n">stud</span><span class="p">,</span><span class="n">t</span><span class="p">),</span> <span class="n">t</span><span class="p">,</span> <span class="mi">0</span><span class="p">,</span> <span class="mi">20</span> <span class="p">)</span></div><div class='line' id='LC14'><span class="k">print</span> <span class="s">&quot;_&quot;</span><span class="o">*</span><span class="mi">50</span></div><div class='line' id='LC15'><span class="k">print</span><span class="p">(</span><span class="n">s</span><span class="p">)</span></div><div class='line' id='LC16'><br/></div></pre></div>
              
            
          </td>
        </tr>
      </table>
    
  </div>


          </div>
        </div>
      </div>
    </div>
  

  </div>


<div class="frame frame-loading" style="display:none;">
  <img src="https://a248.e.akamai.net/assets.github.com/images/modules/ajax/big_spinner_336699.gif" height="32" width="32">
</div>

    </div>
  
      
    </div>

    <div id="footer" class="clearfix">
      <div class="site">
        
          <div class="sponsor">
            <a href="http://www.rackspace.com" class="logo">
              <img alt="Dedicated Server" height="36" src="https://a248.e.akamai.net/assets.github.com/images/modules/footer/rackspace_logo.png?v2" width="38" />
            </a>
            Powered by the <a href="http://www.rackspace.com ">Dedicated
            Servers</a> and<br/> <a href="http://www.rackspacecloud.com">Cloud
            Computing</a> of Rackspace Hosting<span>&reg;</span>
          </div>
        

        <ul class="links">
          
            <li class="blog"><a href="https://github.com/blog">Blog</a></li>
            <li><a href="https://github.com/about">About</a></li>
            <li><a href="https://github.com/contact">Contact &amp; Support</a></li>
            <li><a href="https://github.com/training">Training</a></li>
            <li><a href="http://jobs.github.com">Job Board</a></li>
            <li><a href="http://shop.github.com">Shop</a></li>
            <li><a href="http://developer.github.com">API</a></li>
            <li><a href="http://status.github.com">Status</a></li>
          
        </ul>
        <ul class="sosueme">
          <li class="main">&copy; 2011 <span id="_rrt" title="0.30294s from fe5.rs.github.com">GitHub</span> Inc. All rights reserved.</li>
          <li><a href="/site/terms">Terms of Service</a></li>
          <li><a href="/site/privacy">Privacy</a></li>
          <li><a href="https://github.com/security">Security</a></li>
        </ul>
      </div>
    </div><!-- /#footer -->

    <script>window._auth_token = "a58c21d49dbae0a97a9c153e4e56a698518b585c"</script>
    

<div id="keyboard_shortcuts_pane" class="instapaper_ignore readability-extra" style="display:none">
  <h2>Keyboard Shortcuts <small><a href="#" class="js-see-all-keyboard-shortcuts">(see all)</a></small></h2>

  <div class="columns threecols">
    <div class="column first">
      <h3>Site wide shortcuts</h3>
      <dl class="keyboard-mappings">
        <dt>s</dt>
        <dd>Focus site search</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>?</dt>
        <dd>Bring up this help dialog</dd>
      </dl>
    </div><!-- /.column.first -->

    <div class="column middle" style='display:none'>
      <h3>Commit list</h3>
      <dl class="keyboard-mappings">
        <dt>j</dt>
        <dd>Move selected down</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>k</dt>
        <dd>Move selected up</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>t</dt>
        <dd>Open tree</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>p</dt>
        <dd>Open parent</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>c <em>or</em> o <em>or</em> enter</dt>
        <dd>Open commit</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>y</dt>
        <dd>Expand URL to its canonical form</dd>
      </dl>
    </div><!-- /.column.first -->

    <div class="column last" style='display:none'>
      <h3>Pull request list</h3>
      <dl class="keyboard-mappings">
        <dt>j</dt>
        <dd>Move selected down</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>k</dt>
        <dd>Move selected up</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>o <em>or</em> enter</dt>
        <dd>Open issue</dd>
      </dl>
    </div><!-- /.columns.last -->

  </div><!-- /.columns.equacols -->

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Issues</h3>

    <div class="columns threecols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt>j</dt>
          <dd>Move selected down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>k</dt>
          <dd>Move selected up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>x</dt>
          <dd>Toggle select target</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>o <em>or</em> enter</dt>
          <dd>Open issue</dd>
        </dl>
      </div><!-- /.column.first -->
      <div class="column middle">
        <dl class="keyboard-mappings">
          <dt>I</dt>
          <dd>Mark selected as read</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>U</dt>
          <dd>Mark selected as unread</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>e</dt>
          <dd>Close selected</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>y</dt>
          <dd>Remove selected from view</dd>
        </dl>
      </div><!-- /.column.middle -->
      <div class="column last">
        <dl class="keyboard-mappings">
          <dt>c</dt>
          <dd>Create issue</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>l</dt>
          <dd>Create label</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>i</dt>
          <dd>Back to inbox</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>u</dt>
          <dd>Back to issues</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>/</dt>
          <dd>Focus issues search</dd>
        </dl>
      </div>
    </div>
  </div>

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Network Graph</h3>
    <div class="columns equacols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt><span class="badmono">←</span> <em>or</em> h</dt>
          <dd>Scroll left</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">→</span> <em>or</em> l</dt>
          <dd>Scroll right</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">↑</span> <em>or</em> k</dt>
          <dd>Scroll up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">↓</span> <em>or</em> j</dt>
          <dd>Scroll down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>t</dt>
          <dd>Toggle visibility of head labels</dd>
        </dl>
      </div><!-- /.column.first -->
      <div class="column last">
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">←</span> <em>or</em> shift h</dt>
          <dd>Scroll all the way left</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">→</span> <em>or</em> shift l</dt>
          <dd>Scroll all the way right</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">↑</span> <em>or</em> shift k</dt>
          <dd>Scroll all the way up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">↓</span> <em>or</em> shift j</dt>
          <dd>Scroll all the way down</dd>
        </dl>
      </div><!-- /.column.last -->
    </div>
  </div>

  <div >
    <div class="rule"></div>
    <div class="columns threecols">
      <div class="column first" >
        <h3>Source Code Browsing</h3>
        <dl class="keyboard-mappings">
          <dt>t</dt>
          <dd>Activates the file finder</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>l</dt>
          <dd>Jump to line</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>y</dt>
          <dd>Expand URL to its canonical form</dd>
        </dl>
      </div>
    </div>
  </div>
</div>

    <div id="markdown-help" class="instapaper_ignore readability-extra">
  <h2>Markdown Cheat Sheet</h2>

  <div class="cheatsheet-content">

  <div class="mod">
    <div class="col">
      <h3>Format Text</h3>
      <p>Headers</p>
      <pre>
# This is an &lt;h1&gt; tag
## This is an &lt;h2&gt; tag
###### This is an &lt;h6&gt; tag</pre>
     <p>Text styles</p>
     <pre>
*This text will be italic*
_This will also be italic_
**This text will be bold**
__This will also be bold__

*You **can** combine them*
</pre>
    </div>
    <div class="col">
      <h3>Lists</h3>
      <p>Unordered</p>
      <pre>
* Item 1
* Item 2
  * Item 2a
  * Item 2b</pre>
     <p>Ordered</p>
     <pre>
1. Item 1
2. Item 2
3. Item 3
   * Item 3a
   * Item 3b</pre>
    </div>
    <div class="col">
      <h3>Miscellaneous</h3>
      <p>Images</p>
      <pre>
![GitHub Logo](/images/logo.png)
Format: ![Alt Text](url)
</pre>
     <p>Links</p>
     <pre>
http://github.com - automatic!
[GitHub](http://github.com)</pre>
<p>Blockquotes</p>
     <pre>
As Kanye West said:
> We're living the future so
> the present is our past.
</pre>
    </div>
  </div>
  <div class="rule"></div>

  <h3>Code Examples in Markdown</h3>
  <div class="col">
      <p>Syntax highlighting with <a href="http://github.github.com/github-flavored-markdown/" title="GitHub Flavored Markdown">GFM</a></p>
      <pre>
```javascript
function fancyAlert(arg) {
  if(arg) {
    $.facebox({div:'#foo'})
  }
}
```</pre>
    </div>
    <div class="col">
      <p>Or, indent your code 4 spaces</p>
      <pre>
Here is a Python code example
without syntax highlighting:

    def foo:
      if not bar:
        return true</pre>
    </div>
    <div class="col">
      <p>Inline code for comments</p>
      <pre>
I think you should use an
`&lt;addr&gt;` element here instead.</pre>
    </div>
  </div>

  </div>
</div>
    

    <!--[if IE 8]>
    <script type="text/javascript" charset="utf-8">
      $(document.body).addClass("ie8")
    </script>
    <![endif]-->

    <!--[if IE 7]>
    <script type="text/javascript" charset="utf-8">
      $(document.body).addClass("ie7")
    </script>
    <![endif]-->

    
    
    
    <script type="text/javascript">(function(){var d=document;var e=d.createElement("script");e.async=true;e.src="https://d1ros97qkrwjf5.cloudfront.net/12/eum/rum.js	";e.type="text/javascript";var s=d.getElementsByTagName("script")[0];s.parentNode.insertBefore(e,s);})();NREUMQ.push(["nrf2","beacon-1.newrelic.com","2f94e4d8c2",64799,"dw1bEBZcX1RWRhoEClsAGhcMXEQ=",0.0,299,new Date().getTime()])</script>
  </body>
</html>

