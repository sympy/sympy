document.addEventListener("DOMContentLoaded", function () {

    let toggleButton = document.getElementById("sidebar-toggle");
    if (!toggleButton) {
      toggleButton = document.createElement("button");
      toggleButton.id = "sidebar-toggle";
      toggleButton.setAttribute("aria-label", "Toggle sidebar");
      toggleButton.innerHTML = "☰";
      toggleButton.style.position = "fixed";
      toggleButton.style.top = "1rem";
      toggleButton.style.left = "1rem";
      toggleButton.style.zIndex = "9999";
      toggleButton.style.fontSize = "2rem"
      document.body.insertAdjacentElement("afterbegin", toggleButton);
    }
    
    // make the sidebar provided by furro collapsable
    const sidebar = document.querySelector("aside.sidebar-drawer");
    if (toggleButton && sidebar) {
      toggleButton.addEventListener("click", function () {
        sidebar.classList.toggle("collapsed");
      });
    } else {
      console.error("Toggle button or sidebar-drawer not found");
    }
  });